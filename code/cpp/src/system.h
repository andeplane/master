#pragma once

class Cell;
class Sorter;
class Grid;
class DSMC_IO;
class Random;
class Settings;
class UnitConverter;
class DSMCTimer;
class MoleculeMover;

#include <threadcontrol.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cinifile.h>

using namespace std;

class System {
private:
    void init_positions();
    void init_velocities();
    void init_molecules();
    void init_cells();
	void move();
    void init_randoms();
    void collide();
	void accelerate();
    void maintain_pressure_source();
    void maintain_pressure_drain();
    void maintain_pressure();
    bool remove_molecule_in_pressure_reservoir(bool remove_from_source);
    void find_position_in_reservoirs(double *r, bool find_position_in_source);
    void add_molecule_in_pressure_reservoirs(bool add_in_source);

public:
    DSMC_IO *io;
    DSMCTimer *timer;
    Grid *world_grid;
    Settings *settings;
    UnitConverter * unit_converter;
    Random *rnd;
    MoleculeMover *mover;

    double reservoir_size;
    double grid_origo_x, grid_origo_y, grid_origo_z;
    double length[3];
    double half_length[3];
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double t;
    double temperature;
    double diam, density;
    double wall_temperature;
    double cell_length_x;
    double cell_length_y;
    double cell_length_z;
    double porosity_global;
    double volume;
    int    num_molecules_global;

    unsigned long collisions;
	int steps;
    int myid;
    int cells_x, cells_y, cells_z;
    int num_cells_vector[3];

    ThreadControl thread_control;

    void initialize(Settings *settings_, int myid_);
	void step();
    System() { }
};
