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

void Molecule::move(double dt, Random *rnd) {

    double L = system->L;

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

    r(0) = fmod(r(0)+10*L,L);
    r(1) = fmod(r(1)+10*L,L);
}
