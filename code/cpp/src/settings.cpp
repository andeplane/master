#include "settings.h"

Settings::Settings(string filename) {
    ini_file.load(filename);

    load_previous_state = ini_file.getbool("load_previous_state");
    create_movie = ini_file.getbool("create_movie");
    statistics_interval = ini_file.getint("statistics_interval");
    movie_every_n_frame = ini_file.getint("movie_every_n_frame");
    number_of_molecules = ini_file.getint("N");
    timesteps = ini_file.getint("timesteps");

    cells_per_node_x = ini_file.getint("cells_per_node_x");
    cells_per_node_y = ini_file.getint("cells_per_node_y");
    cells_per_node_z = ini_file.getint("cells_per_node_z");
    nodes_x = ini_file.getint("nodes_x");
    nodes_y = ini_file.getint("nodes_y");
    nodes_z = ini_file.getint("nodes_z");
    gravity_direction = ini_file.getint("gravity_direction");

    viscosity = ini_file.getdouble("viscosity");
    mass = ini_file.getdouble("mass");
    temperature = ini_file.getdouble("temperature");
    wall_temperature = ini_file.getdouble("wall_temperature");
    gravity = ini_file.getdouble("gravity");
    density = ini_file.getdouble("density");
    diam = ini_file.getdouble("diam");
    Lx = ini_file.getdouble("Lx");
    Ly = ini_file.getdouble("Ly");
    Lz = ini_file.getdouble("Lz");
    L_reservoir_x = ini_file.getdouble("L_reservoir_x");
    L_reservoir_y = ini_file.getdouble("L_reservoir_y");
    L_reservoir_z = ini_file.getdouble("L_reservoir_z");

    dt = ini_file.getdouble("dt");
}
