#include "settings.h"

Settings::Settings(string filename) {
    ini_file.load(filename);

    load_previous_state = ini_file.getbool("load_previous_state");
    create_movie = ini_file.getbool("create_movie");
    statistics_interval = ini_file.getint("statistics_interval");
    movie_every_n_frame = ini_file.getint("movie_every_n_frame");
    number_of_particles = ini_file.getint("N");
    timesteps = ini_file.getint("timesteps");

    cells_per_node_x = ini_file.getint("cells_per_node_x");
    cells_per_node_y = ini_file.getint("cells_per_node_y");
    cells_per_node_z = ini_file.getint("cells_per_node_z");
    nodes_x = ini_file.getint("nodes_x");
    nodes_y = ini_file.getint("nodes_y");
    nodes_z = ini_file.getint("nodes_z");

    temperature = ini_file.getdouble("temperature");
    wall_temperature = ini_file.getdouble("wall_temperature");
    acceleration = ini_file.getdouble("acceleration");
    max_x_acceleration = ini_file.getdouble("max_x_acceleration");
    density = ini_file.getdouble("density");
    diam = ini_file.getdouble("diam");
    Lx = ini_file.getdouble("Lx");
    Ly = ini_file.getdouble("Ly");
    Lz = ini_file.getdouble("Lz");
    dt = ini_file.getdouble("dt");
    cells_x = ini_file.getint("cells_x");
    cells_y = ini_file.getint("cells_y");
    cells_z = ini_file.getint("cells_z");
}
