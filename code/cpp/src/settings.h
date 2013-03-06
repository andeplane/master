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
    int cells_per_node_x;
    int cells_per_node_y;
    int cells_per_node_z;
    int nodes_x;
    int nodes_y;
    int nodes_z;
    double temperature;
    double wall_temperature;
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
