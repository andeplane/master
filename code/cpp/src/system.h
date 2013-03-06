#pragma once

class Molecule;
class Cell;
class Sorter;
class Grid;
class DSMC_IO;
class Random;
class Settings;
class UnitConverter;

#include <iostream>
#include <fstream>
#include <vector>
#include <CIniFile.h>
#include <settings.h>


using namespace std;

class System {
private:
    void init_positions();
    void init_velocities();
    void init_molecules();
    void init_cells();
	void move();
    void init_randoms();
	int  collide();
	void accelerate();
public:
    vector<Molecule*>molecules;
    vector< vector< vector<Cell*> > > cells;

    DSMC_IO *io;
    Grid *world_grid;
    Settings *settings;
    UnitConverter * unit_converter;

    Random *rnd;

	int N; 			// Number of molecules

    double Lx;
    double Ly;
    double Lz;
    double acceleration;
    double max_x_acceleration;
	double volume;
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double t;
    double temperature;
    double mass, diam, density;
    double wall_temperature;
	double *time_consumption;

    double *positions;
    double *velocities;
    double *initial_positions;

	int collisions;
	int steps;

	Sorter *sorter;

    void initialize(Settings *settings_);
	void step();
    System() { }
};
