#include "Molecule.h"
#include <iostream>
#include <math.h>
#include "lib.h"
#include "Wall.h"
#include "omp.h"
using namespace std;

Molecule::Molecule(System *system) {
	this->r = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);

	this->atoms = 1;
	this->type = 0;
	this->system = system;
}

static int fuck = 0;

inline int sign(double a) {
	return a >= 0 ? 1 : -1;
}

void Molecule::move(double dt) {
	double L = this->system->L;
	
	vec v = this->v;
	double y_old = this->r(1);
	this->r += v*dt;

	int outsideWallIndex = this->system->walls[0]->isMoleculeOutside(this) + 2*this->system->walls[1]->isMoleculeOutside(this);
	double tau = 1;

	if(outsideWallIndex) {
		long *idum = this->system->idums[omp_get_thread_num()];

		Wall *wall = this->system->walls[--outsideWallIndex]; // Decrease from [1,2] to [0,1] first
		int direction = wall->upper ? -1 : 1;

		// We are outside the box
		tau = wall->timeUntilCollision(y_old,v(1));
		
		r += v*(tau-dt);
		
		v(0) = sqrt(3.0/2*wall->T)*system->rand_gauss(idum) + wall->v_x; // Add wall velocity
      	v(1) = direction*sqrt(-3.0/2*wall->T*log(ran0(idum)));

      	this->r += v*(dt-tau);
	}

	this->v = v;
	this->r(0) = fmod(this->r(0)+10*L,L);
}