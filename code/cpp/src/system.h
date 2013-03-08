#pragma once

class Molecule;
class Cell;
class Sorter;
class Grid;
class DSMC_IO;
class Random;
class Settings;
class UnitConverter;

#include <threadcontrol.h>
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
    vector< vector< vector<Cell*> > > cells;

    DSMC_IO *io;
    Grid *world_grid;
    Settings *settings;
    UnitConverter * unit_converter;

    Random *rnd;
    double Lx;
    double Ly;
    double Lz;
    double acceleration;
    double max_x_acceleration;
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double t;
    double temperature;
    double mass, diam, density;
    double wall_temperature;
    double cell_length_x;
    double cell_length_y;
    double cell_length_z;
    double porosity_global;
    double volume;
    int    num_particles_global;

	int collisions;
	int steps;
    int myid;
    int cells_x, cells_y, cells_z;

    ThreadControl thread_control;

    void initialize(Settings *settings_, int myid_);
	void step();
    System() { }
};
