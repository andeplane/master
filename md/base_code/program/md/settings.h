#pragma once

#include <cinifile.h>

class Settings
{
public:
    Settings(string filename);
    CIniFile ini_file;
    double temperature;
    double dt;
    double FCC_b;
    double r_cut;
    double thermostat_relaxation_time;
    double gravity_force;
    double mass;
    int gravity_direction;
    int unit_cells_x;
    int unit_cells_y;
    int unit_cells_z;
    int nodes_x;
    int nodes_y;
    int nodes_z;
    int timesteps;
    int movie_every_n_frame;
    int statistics_interval;
    bool create_movie;
    bool load_state;
    bool thermostat_enabled;
    bool many_frozen_atoms;

};
