#pragma once

#include <cinifile.h>

class Settings
{
public:
    Settings(string filename);
    CIniFile ini_file;
    bool load_previous_state;
    bool create_movie;
    bool maintain_pressure;
    int statistics_interval;
    int movie_every_n_frame;
    int atoms_per_molecule;
    int movie_molecules;
    int timesteps;
    int cells_x;
    int cells_y;
    int cells_z;
    int gravity_direction;
    double viscosity;
    double mass;
    double temperature;
    double wall_temperature;
    double gravity;
    double density;
    double diam;
    double Lx;
    double Ly;
    double Lz;
    double reservoir_fraction;
    double pressure_source;
    double pressure_drain;
    double dt;
};
