#ifndef SYSTEM_H
#define SYSTEM_H

class Molecule;
class Cell;
class Sorter;

#include <fstream>
#include "Molecule.h"
#include "Cell.h"
#include "Sorter.h"

using namespace std;
using namespace arma;

class System {
private:
	void initialize();
	void initPositions();
	void initVelocities();
	void initMolecules();
	void initCells();
	void move();
	int  collide();
public:
	Molecule **molecules;
	Cell **cells;
	int N; 			// Number of molecules
	double L;
	double volume;
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double vwall;
	double tau;
	double coeff;
	double t;
	double T;
	int numberOfCells;
	int cellsPerDimension;
	int collisions;
	int steps;

	vec wallStrikes;
	vec delta_v;
	
	long *idum;

	Sorter *sorter;

	void printPositionsToFile(FILE *file);

	void step();
	System(int N=500, double T=300.0);

};


#endif