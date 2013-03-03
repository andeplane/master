#include "defines.h"
#include "molecule.h"
#include <iostream>
#include <math.h>
#include <CVector.h>
using namespace std;

Molecule::Molecule(System *_system) {
    r = zeros<vec> (3,1);
    initial_r = r;
    v = zeros<vec> (3,1);
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
    temp_vector = initial_r-r;
    return dot(temp_vector,temp_vector);
}

inline void Molecule::do_move(const double &dt) {
    r(0) += v(0)*dt;
    r(1) += v(1)*dt;
    r(2) += v(2)*dt;
    fixR();
}

inline void Molecule::fixR() {
    if(r(0) > system->width)  { r(0) -= system->width; initial_r(0) -= system->width; }
    else if(r(0) < 0)         { r(0) += system->width; initial_r(0) += system->width; }

    if(r(1) > system->height) { r(1) -= system->height; initial_r(1) -= system->height; }
    else if(r(1) < 0)         { r(1) += system->height; initial_r(1) += system->height; }
}

void Molecule::move(double dt, Random *rnd, int depth) {
    if(!active) return;

    double tau = dt;
    do_move(dt);

    GridPoint *point = system->world_grid->get_grid_point(r);

    // We have to calculate time until collision
    if(point->is_wall) {
        if(!point->is_wall_boundary) {
            int count = 0;
            while(true) {
                point = system->world_grid->get_grid_point(r);
                if(point->is_wall_boundary) {
                    dt -= tau;

                    while(system->world_grid->get_grid_point(r)->is_wall) {
                        dt += 0.1*tau;
                        do_move(-0.1*tau);
                    }
                    break;
                }

                if(point->is_wall) {
                    do_move(-tau);
                    tau /= 2;
                }
                else {
                    dt -= tau;
                }

                if(++count > 100) {
                    active = false;
                    cout << "Trouble with molecule " << index << " at (i,j)=(" << point->i << "," << point->j << "). Check your world." << endl;
                    return;
                }
                do_move(tau);
            }
        }
        else {
            // This was actually the boundary2
            // cout << "Particle " << index << " is at the boundary!" << endl;
            int count = 0;
            while(system->world_grid->get_grid_point(r)->is_wall) {
                dt += 0.1*tau;
                do_move(-0.1*tau);
                if(++count > 100) {
                    active = false;
                    cout << "Another kind of trouble with molecule " << index << " at (i,j)=(" << point->i << "," << point->j << "). Check your world." << endl;
                    return;
                }
            }
        }

        double v_normal = sqrt(-6.0/2*point->T*log(rnd->nextDouble()));
        double v_tangent = sqrt(3.0/2*point->T)*rnd->nextGauss();

        v(0) = v_normal*point->normal.x + v_tangent*point->tangent.x;
        v(1) = v_normal*point->normal.y + v_tangent*point->tangent.y;

    }
    else dt = 0;

    if(dt > 1e-10 && depth < 5) {
        move(dt,rnd,depth+1);
    }
}

void Molecule::collide_with(Molecule *m, Random *rnd, const double &cr) {
    double cos_th, sin_th;
    vcm(0) = 0.5*(v(0)+m->v(0));
    vcm(1) = 0.5*(v(1)+m->v(1));
    vcm(2) = 0.5*(v(2)+m->v(2));

    cos_th = 1.0 - 2.0*rnd->nextDouble();      // Cosine and sine of
    sin_th = sqrt(1.0 - cos_th*cos_th); // collision angle theta

    vrel(0) = cr*cos_th;             // Compute post-collision
    vrel(1) = cr*sin_th;             // relative velocity

    temp_vector(0) = atoms*(vcm(0) + 0.5*vrel(0)-v(0));
    temp_vector(1) = atoms*(vcm(1) + 0.5*vrel(1)-v(1));
    temp_vector(2) = atoms*(vcm(2) + 0.5*vrel(2)-v(2));

    v(0) = vcm(0) + 0.5*vrel(0);
    v(1) = vcm(1) + 0.5*vrel(1);
    v(2) = vcm(2) + 0.5*vrel(2);

    m->v(0) = vcm(0) - 0.5*vrel(0);
    m->v(1) = vcm(1) - 0.5*vrel(1);
    m->v(2) = vcm(2) - 0.5*vrel(2);
}
