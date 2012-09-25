#include "Molecule.h"
#include <iostream>
#include <math.h>

using namespace std;

Molecule::Molecule(System *system) {
	this->r = zeros<vec> (3,1);
	this->r_old = zeros<vec> (3,1);
	this->v = zeros<vec> (3,1);
	

	this->mass = 1; // 39.948;         // MD units
	this->type = 0;
	this->system = system;
}

void Molecule::move(double tau) {
	double L = this->system->L;

	this->r_old = this->r;
	this->r += this->v*tau;
	this->r(1) = fmod(this->r(1)+10*L,L);
	this->r(2) = fmod(this->r(2)+10*L,L);
}