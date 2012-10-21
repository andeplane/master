#ifndef SYSTEM_H
#define SYSTEM_H

class Molecule;
class Cell;
class Sorter;
class Wall;

#include <fstream>
#include "Molecule.h"
#include "Cell.h"
#include "Sorter.h"
#include "Wall.h"
#include "Random.h"
#include <vector>
#include <CIniFile.h>

using namespace std;
using namespace arma;

class System {
private:
	void initPositions();
	void initVelocities();
	void initMolecules();
	void initCells();
	void initWalls();
	void move();
	int  collide();
	void accelerate();
public:
    Molecule **molecules;
    Cell ***cells;

	Wall **walls;
    Random **randoms;

	int N; 			// Number of molecules

    double width;
    double height;

	double volume;
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double coeff;
	double t;
	double T;
	double *time_consumption;
	int numberOfCells;
    int threads;
	int collisions;
	int steps;

    int cells_x;
    int cells_y;

	Sorter *sorter;

	void printPositionsToFile(FILE *file);
    void initialize(CIniFile &ini);
	void step();
	double rand_gauss(long *idum);
    System() { }
};


#endif
