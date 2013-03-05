#include "defines.h"
#include "molecule.h"
#include <iostream>
#include <math.h>
#include <CVector.h>
using namespace std;

Molecule::Molecule(System *_system) {
    vrel = zeros<vec>(3,1);
    vcm = zeros<vec>(3,1);
    temp_vector = zeros<vec>(3,1);

    mass = 1;
    atoms = 1;
    type = 0;
    system = _system;
    type = "H";
    information_carrier = 0;
    active = true;
}

inline int sign(double a) {
    return a >= 0 ? 1 : -1;
}

double Molecule::squared_distance_from_initial_position() {
    return pow(initial_r[0] - r[0],2)+pow(initial_r[1] - r[1],2)+pow(initial_r[2] - r[2],2);
}

inline void Molecule::do_move(const double &dt) {
    r[0] += v[0]*dt;
    r[1] += v[1]*dt;
    r[2] += v[2]*dt;
    fixR();
}

inline void Molecule::fixR() {
    if(r[0] > system->Lx)  { r[0] -= system->Lx; initial_r[0] -= system->Lx; }
    else if(r[0] < 0)         { r[0] += system->Lx; initial_r[0] += system->Lx; }

    if(r[1] > system->Ly) { r[1] -= system->Ly; initial_r[1] -= system->Ly; }
    else if(r[1] < 0)         { r[1] += system->Ly; initial_r[1] += system->Ly; }

    if(r[2] > system->Lz) { r[2] -= system->Lz; initial_r[2] -= system->Lz; }
    else if(r[2] < 0)         { r[2] += system->Lz; initial_r[2] += system->Lz; }
}

void Molecule::move(double dt, Random *rnd, int depth) {
    if(!active) return;
    // cout << "################# I AM MOLECULE at " << endl;
    // cout << r[0] << " " << r[1] << " " << r[2] << endl;

    double tau = dt;
    do_move(dt);

    unsigned char *voxels = system->world_grid->voxels;
    int idx = system->world_grid->get_index_of_voxel(r);

    // We have to calculate time until collision
    if(voxels[idx]>=voxel_type_wall) { // Is wall
        if(voxels[idx]!=voxel_type_boundary) { // Not boundary
            int count = 0;
            while(true) {
                idx = system->world_grid->get_index_of_voxel(r);
                if(voxels[idx]==voxel_type_boundary) {
                    // cout << r[0] << " " << r[1] << " " << r[2] << endl;
                    // cout << "Voxel is boundary, accepted this move, removing tau from dt" << endl;
                    dt -= tau;
                    while(system->world_grid->get_voxel(r)[0]>=voxel_type_wall) {
                        // cout << "Voxel is still boundary, going a bit back..." << endl;
                        dt += 0.1*tau;
                        do_move(-0.1*tau);
                    }
                    break;
                }

                if(voxels[idx]>=voxel_type_wall) {
                    // cout << r[0] << " " << r[1] << " " << r[2] << endl;
                    // cout << "Voxel is wall, going back..." << endl;
                    do_move(-tau);
                    tau /= 2;
                }
                else {
                    // cout << r[0] << " " << r[1] << " " << r[2] << endl;
                    // cout << "Accepted this move, removing tau from dt" << endl;
                    dt -= tau;
                }

                if(++count > 100) {
                    active = false;
                    // cout << "I have trouble..." << endl;
                    // cout << "Trouble with molecule " << index << " at (i,j)=(" << point->i << "," << point->j << "). Check your world." << endl;
                    exit(0);
                }

                do_move(tau);
            }
        }
        else {
            int count = 0;
            while(*system->world_grid->get_voxel(r)>=voxel_type_wall) {
                dt += 0.1*tau;
                do_move(-0.1*tau);
                if(++count > 100) {
                    active = false;
                    // cout << "I have another kind of trouble" << endl;
                    // cout << "Another kind of trouble with molecule " << index << " at (i,j)=(" << point->i << "," << point->j << "). Check your world." << endl;
                    return;
                }
            }
        }

        double v_normal   = sqrt(-2*system->wall_temperature*log(rnd->nextDouble()));
        double v_tangent1 = sqrt(system->wall_temperature)*rnd->nextGauss();
        double v_tangent2 = sqrt(system->wall_temperature)*rnd->nextGauss();

        // Normal vector
        double n_x = system->world_grid->normals[3*idx + 0];
        double n_y = system->world_grid->normals[3*idx + 1];
        double n_z = system->world_grid->normals[3*idx + 2];

        // Tangent vector 1
        double t1_x = system->world_grid->tangents1[3*idx + 0];
        double t1_y = system->world_grid->tangents1[3*idx + 1];
        double t1_z = system->world_grid->tangents1[3*idx + 2];

        // Tangent vector 2
        double t2_x = system->world_grid->tangents2[3*idx + 0];
        double t2_y = system->world_grid->tangents2[3*idx + 1];
        double t2_z = system->world_grid->tangents2[3*idx + 2];

        v[0] = v_normal*n_x + v_tangent1*t1_x + v_tangent2*t2_x;
        v[1] = v_normal*n_y + v_tangent1*t1_y + v_tangent2*t2_y;
        v[2] = v_normal*n_z + v_tangent1*t1_z + v_tangent2*t2_z;
    }
    else dt = 0;

    if(dt > 1e-10 && depth < 5) {
        move(dt,rnd,depth+1);
    }
}

void Molecule::collide_with(Molecule *m, Random *rnd, const double &cr) {
    double cos_th, sin_th, phi, vcmx, vcmy, vcmz, vrelx, vrely, vrelz, tmpx, tmpy, tmpz;
    vcmx  = 0.5*(v[0] + m->v[0]);
    vcmy  = 0.5*(v[1] + m->v[1]);
    vcmz  = 0.5*(v[2] + m->v[2]);

    cos_th = 1.0 - 2.0*rnd->nextDouble();      // Cosine and sine of
    sin_th = sqrt(1.0 - cos_th*cos_th); // collision angle theta
    phi = 2*M_PI*rnd->nextDouble();


    vrelx = cr*cos_th;                   // Compute post-collision
    vrely = cr*sin_th*cos(phi);          // relative velocity
    vrelz = cr*sin_th*sin(phi);

    tmpx = atoms*(vcmx + 0.5*vrelx-v[0]);
    tmpy = atoms*(vcmy + 0.5*vrely-v[1]);
    tmpz = atoms*(vcmz + 0.5*vrelz-v[2]);

    v[0] = vcmx + 0.5*vrelx;
    v[1] = vcmy + 0.5*vrely;
    v[2] = vcmz + 0.5*vrelz;

    m->v[0] = vcmx - 0.5*vrelx;
    m->v[1] = vcmy - 0.5*vrely;
    m->v[2] = vcmz - 0.5*vrelz;
}
