#include "Molecule.h"
#include <iostream>
#include <math.h>
#include "Wall.h"
#include "omp.h"
using namespace std;

Molecule::Molecule(System *_system) {
    r = zeros<vec> (3,1);
    v = zeros<vec> (3,1);

    atoms = 1;
    type = 0;
    system = _system;
}

inline int sign(double a) {
	return a >= 0 ? 1 : -1;
}

void Molecule::move(double dt, Random *rnd, int depth) {
    move_old(dt,rnd);
    return;

    int resolution = 5;

    double tau = dt;
    r += v*dt;

    bool didCollide = system->world_grid->get_grid_point(r(0),r(1));

    // We have to calculate time until collision
    if(didCollide) {
        for(int d=1;d<resolution;d++) {
            // Go back in time, try with another timestep
            r -= v*tau;
            tau = d*dt/resolution;
            r += v*tau;

            didCollide = system->world_grid->get_grid_point(r(0),r(1));
            if(didCollide) {
                r -= v*tau;
                tau = (d-1)*dt/resolution;
                r += v*tau;

                break;
            }
        }
        /*
        int sign = 2*(v(0)*v(1)) - 1;
        double theta = -3*sign*M_PI/2;
        v(0) = v(0)*cos(theta) - v(1)*sin(theta);
        v(1) = v(0)*sin(theta) + v(1)*cos(theta);
        */
        v = -v;
    }
    if(dt-tau > 1e-8)
        move(dt-tau,rnd,depth+1);

    if(depth == 0) {
        r(0) = fmod(r(0)+10*system->width,system->width);
        r(1) = fmod(r(1)+10*system->height,system->height);
    }
}

void Molecule::move_old(double dt, Random *rnd) {
    double y_old = r(1);
    r += v*dt;

    int outsideWallIndex = system->walls[0]->isMoleculeOutside(this) + 2*system->walls[1]->isMoleculeOutside(this);
	double tau = 1;

	if(outsideWallIndex) {
        Wall *wall = system->walls[--outsideWallIndex]; // Decrease from [1,2] to [0,1] first
		int direction = wall->upper ? -1 : 1;

		// We are outside the box
		tau = wall->timeUntilCollision(y_old,v(1));
		
        r += v*(tau-dt);

        v(0) = sqrt(3.0/2*wall->T)*rnd->nextGauss() + wall->v_x;
        v(1) = direction*sqrt(-6.0/2*wall->T*log(rnd->nextDouble()));

        r += v*(dt-tau);
	}

    r(0) = fmod(r(0)+10*system->width,system->width);
    r(1) = fmod(r(1)+10*system->height,system->height);
}
