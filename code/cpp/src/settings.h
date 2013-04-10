#pragma once

#include <cinifile.h>

class Settings
{
public:
    Settings(string filename);
    CIniFile ini_file;
    bool load_previous_state;
    bool create_movie;
    int statistics_interval;
    int movie_every_n_frame;
    int number_of_molecules;
    int timesteps;
    int cells_per_node_x;
    int cells_per_node_y;
    int cells_per_node_z;
    int nodes_x;
    int nodes_y;
    int nodes_z;
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
    double L_reservoir_x;
    double L_reservoir_y;
    double L_reservoir_z;
    double dt;
};
