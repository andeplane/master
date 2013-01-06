#include "defines.h"
#include "molecule.h"
#include <iostream>
#include <math.h>
#include "omp.h"
#include <CVector.h>
using namespace std;

Molecule::Molecule(System *_system) {
    r = zeros<vec> (3,1);
    v = zeros<vec> (3,1);

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

inline void Molecule::addR(vec dr) {
    r += dr;
    fixR();
}

inline void Molecule::fixR() {
    if(r(0) > system->width) r(0) -= system->width;
    else if(r(0) < 0)        r(0) += system->width;
    if(r(1) > system->height) r(1) -= system->height;
    else if(r(1) < 0)        r(1) += system->height;
}

#ifdef VISCOSITY
    void Molecule::move(double dt, Random *rnd, int depth) {
        r += v*dt;

        if(r(1) < 0) {
            r -= v*dt;
            double tau = r(1)/abs(v(1));
            r += v*tau;

            v(1) = sqrt(-6.0/2*system->wall_temperature*log(rnd->nextDouble()));
            v(0) = sqrt(3.0/2*system->wall_temperature)*rnd->nextGauss();

            r += v*(dt - tau);
        }
        else if(r(1) > system->height) {
            r -= v*dt;
            double tau = (system->height - r(1))/abs(v(1));
            r += v*tau;

            v(1) = -sqrt(-6.0/2*system->wall_temperature*log(rnd->nextDouble()));
            v(0) = sqrt(3.0/2*system->wall_temperature)*rnd->nextGauss() + VWALL;  // Wall velocity
            r += v*(dt - tau);
        }
        if(r(0) < 0 || r(0) > system->width) {
            r(0) = fmod(r(0) + 10*system->width,system->width);
        }
    }

#else


void Molecule::move(double dt, Random *rnd, int depth) {
    if(!active) return;

    double tau = dt;

    addR(v*dt);

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
                        addR(-0.1*v*tau);
                    }
                    break;
                }

                if(point->is_wall) {
                    addR(-v*tau);
                    tau /= 2;
                }
                else {
                    dt -= tau;
                }

                if(++count > 100) {
                    active = false;
                    r.zeros();
                    v.zeros();
                    cout << "Trouble with molecule " << index << " at (i,j)=(" << point->i << "," << point->j << "). Check your world." << endl;
                    return;
                }

                addR(v*tau);
            }
        }
        else {
            // This was actually the boundary2
            // cout << "Particle " << index << " is at the boundary!" << endl;
            int count = 0;
            while(system->world_grid->get_grid_point(r)->is_wall) {
                dt += 0.1*tau;
                addR(-0.1*v*tau);
                if(++count > 100) {
                    active = false;
                    r.zeros();
                    cout << "Another kind of trouble with molecule " << index << " at (i,j)=(" << point->i << "," << point->j << "). Check your world." << endl;
                    return;
                }
            }
        }

        // double v_normal = sqrt(-6.0/2*point->T*log(rnd->nextDouble()));
        double v_normal = sqrt(-6.0/2*point->T*log(rnd->nextDouble()));
        double v_tangent = sqrt(3.0/2*point->T)*rnd->nextGauss();

        v(0) = v_normal*point->normal.x + v_tangent*point->tangent.x;
        v(1) = v_normal*point->normal.y + v_tangent*point->tangent.y;

    }
    else dt = 0;

    if(dt > 1e-10 && depth < 5)
        move(dt,rnd,depth+1);
}
#endif
