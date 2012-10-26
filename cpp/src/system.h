#pragma once

class Molecule;
class Cell;
class Sorter;
class Wall;
class Grid;

#include <fstream>
#include "molecule.h"
#include "cell.h"
#include "sorter.h"
#include "random.h"
#include <vector>
#include <CIniFile.h>
#include <grid.h>

using namespace std;
using namespace arma;

class System {
private:
	void initPositions();
	void initVelocities();
	void initMolecules();
	void initCells();
	void move();
	int  collide();
	void accelerate();
public:
    Molecule **molecules;
    Cell ***cells;
    Grid *world_grid;
    Grid *initial_world_grid;

    Random **randoms;

	int N; 			// Number of molecules

    double width;
    double height;
    double acceleration;
    double max_x_acceleration;
	double volume;
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double coeff;
	double t;
	double T;
    double wall_temperature;
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
