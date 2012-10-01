#ifndef WALL_H
#define WALL_H

#include "Molecule.h"
#include "System.h"
#include <armadillo>

using namespace arma;
using namespace std;

class Wall {
private:
	
public:
	double y;
	double T;
	double v_x;
	bool   upper;

	Wall(double y, double T, bool upper, double v_x=0.0);
	bool isMoleculeOutside(Molecule *molecule);
	double timeUntilCollision(double y_old, double v_y);
};

#endif