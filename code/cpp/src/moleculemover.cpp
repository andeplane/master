#include <moleculemover.h>
#include <system.h>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <grid.h>
#include <threadcontrol.h>
#include <settings.h>
MoleculeMover::MoleculeMover()
{

}

double sqrt_wall_temp_over_mass = 0;

void MoleculeMover::initialize(System *system_) {
    system = system_;
    voxels = system->world_grid->voxels;
    grid = system->world_grid;
    thread_control = &system->thread_control;
    sqrt_wall_temp_over_mass = sqrt(system->wall_temperature/system->settings->mass);
    count_periodic_x = 0;
    count_periodic_y = 0;
    count_periodic_z = 0;
}

void MoleculeMover::move_molecules(double dt, Random *rnd) {
    for(int n=0;n<thread_control->num_molecules;n++) {
        move_molecule(n,dt,rnd,0);
    }
}

void MoleculeMover::do_move(double *r, double *v, double *r0, const double &dt) {
    r[0] += v[0]*dt;
    r[1] += v[1]*dt;
    r[2] += v[2]*dt;

    if(r[0] > system->Lx)  { r[0] -= system->Lx; r0[0] -= system->Lx; count_periodic_x++; }
    else if(r[0] < 0)         { r[0] += system->Lx; r0[0] += system->Lx; count_periodic_x--; }

    if(r[1] > system->Ly) { r[1] -= system->Ly; r0[1] -= system->Ly; count_periodic_y++; }
    else if(r[1] < 0)         { r[1] += system->Ly; r0[1] += system->Ly; count_periodic_y--; }

    if(r[2] > system->Lz) { r[2] -= system->Lz; r0[2] -= system->Lz; count_periodic_z++; }
    else if(r[2] < 0)         { r[2] += system->Lz; r0[2] += system->Lz; count_periodic_z--; }
}

inline int get_index_of_voxel(double *r, const double &nx_div_lx,const double &ny_div_ly,const double &nz_div_lz, const int &Nx, const int &NyNx) {
    int i =  r[0]*nx_div_lx;
    int j =  r[1]*ny_div_ly;
    int k =  r[2]*nz_div_lz;

    return i + j*Nx + k*NyNx;
}

void MoleculeMover::move_molecule(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;

    double *r = &thread_control->r[3*molecule_index];
    double *r0 = &thread_control->r0[3*molecule_index];
    double *v = &thread_control->v[3*molecule_index];

    do_move(r,v,r0,dt);

    double nx_div_lx = grid->Nx/system->Lx;
    double ny_div_ly = grid->Ny/system->Ly;
    double nz_div_lz = grid->Nz/system->Lz;
    int nx = grid->Nx;
    int nynx = grid->Nx*grid->Ny;

    int idx = get_index_of_voxel(r,nx_div_lx,ny_div_ly,nz_div_lz,nx,nynx);

    // We have to calculate time until collision
    if(voxels[idx]>=voxel_type_wall) { // Is wall
        if(voxels[idx]!=voxel_type_boundary) { // Not boundary
            int count = 0;
            while(true) {
                idx = get_index_of_voxel(r,nx_div_lx,ny_div_ly,nz_div_lz,nx,nynx);
                if(voxels[idx]==voxel_type_boundary) {
                    dt -= tau;
                    while(*grid->get_voxel(r)>=voxel_type_wall) {
                        dt += 0.1*tau;
                        do_move(r,v,r0,-0.1*tau);
                    }
                    break;
                }

                if(voxels[idx]>=voxel_type_wall) {
                    do_move(r,v,r0,-tau);
                    tau /= 2;
                }
                else {
                    dt -= tau;
                }

                if(++count > 100) {
                    cout << "TROUBLE 1" << endl;
                    exit(0);
                }

                do_move(r,v,r0,tau);
            }
        }
        else {
            int count = 0;
            while(*grid->get_voxel(r)>=voxel_type_wall) {
                dt += 0.1*tau;
                do_move(r,v,r0,-0.1*tau);
                if(++count > 100) {
                    cout << "TROUBLE 2" << endl;
                    exit(0);
                }
            }
        }

        double v_normal   = sqrt_wall_temp_over_mass*sqrt(-2*log(rnd->nextDouble()));
        double v_tangent1 = sqrt_wall_temp_over_mass*rnd->nextGauss();
        double v_tangent2 = sqrt_wall_temp_over_mass*rnd->nextGauss();

        // Normal vector
        float n_x = grid->normals[3*idx + 0];
        float n_y = grid->normals[3*idx + 1];
        float n_z = grid->normals[3*idx + 2];

        // Tangent vector 1
        float t1_x = grid->tangents1[3*idx + 0];
        float t1_y = grid->tangents1[3*idx + 1];
        float t1_z = grid->tangents1[3*idx + 2];

        // Tangent vector 2
        float t2_x = grid->tangents2[3*idx + 0];
        float t2_y = grid->tangents2[3*idx + 1];
        float t2_z = grid->tangents2[3*idx + 2];

        v[0] = v_normal*n_x + v_tangent1*t1_x + v_tangent2*t2_x;
        v[1] = v_normal*n_y + v_tangent1*t1_y + v_tangent2*t2_y;
        v[2] = v_normal*n_z + v_tangent1*t1_z + v_tangent2*t2_z;
    }
    else dt = 0;

    if(dt > 1e-5 && depth < 10) {
        move_molecule(molecule_index,dt,rnd,depth+1);
    }
}
