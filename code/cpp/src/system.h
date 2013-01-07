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
    void init_randoms();
	int  collide();
	void accelerate();
public:
    Molecule **molecules;
    Cell ***cells;
    Cell ***load_balanced_cell_list;
    int  *cells_in_list;

    Grid *world_grid;
    Grid *initial_world_grid;

    Random **randoms;

	int N; 			// Number of molecules

    double dvx0;
    double dvx1;

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
    double mass, diam, density;
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
    void read_ini_file(CIniFile &ini);
	double rand_gauss(long *idum);
    System() { }
};
