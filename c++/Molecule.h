#ifndef MOLECULE_H
#define MOLECULE_H
class System;

#include <iostream>
#include <armadillo>
#include "System.h"

using namespace arma;
using namespace std;

class Molecule {
public:
	char   type;
	double mass;

	vec r;
	vec r_old;
	vec v;

	System *system;
	
	Molecule(System *system);
	void addR(vec dr);
	void move(double tau);
};

#endif