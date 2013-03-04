#pragma once

#include <CIniFile.h>

class Settings
{
public:
    Settings(string filename);
    CIniFile ini_file;
    bool load_previous_state;
    bool create_movie;
    int statistics_interval;
    int movie_every_n_frame;
    int number_of_particles;
    int timesteps;
    double temperature;
    double wall_temperature;
    double dt_factor;
    double acceleration;
    double max_x_acceleration;
    double density;
    double diam;
    double Lx;
    double Ly;
    double Lz;
    double dt;
    int cells_x;
    int cells_y;
    int cells_z;

};
