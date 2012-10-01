#include "Molecule.h"
#include <iostream>
#include <math.h>
using namespace std;

Molecule::Molecule(System *system) {
	this->r = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);

	this->atoms = 1;
	this->type = 0;
	this->system = system;
}

void Molecule::move(double dt) {
	double L = this->system->L;
	this->r += this->v*dt;
	
	this->r(0) = fmod(this->r(0)+10*L,L);
	this->r(1) = fmod(this->r(1)+10*L,L);
}