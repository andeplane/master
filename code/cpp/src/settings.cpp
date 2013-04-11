#include "settings.h"

Settings::Settings(string filename) {
    ini_file.load(filename);

    load_previous_state = ini_file.getbool("load_previous_state");
    create_movie = ini_file.getbool("create_movie");
    maintain_pressure = ini_file.getbool("maintain_pressure");
    statistics_interval = ini_file.getint("statistics_interval");
    movie_every_n_frame = ini_file.getint("movie_every_n_frame");
    number_of_molecules = ini_file.getint("N");
    timesteps = ini_file.getint("timesteps");

    cells_x = ini_file.getint("cells_x");
    cells_y = ini_file.getint("cells_y");
    cells_z = ini_file.getint("cells_z");
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
    reservoir_fraction = ini_file.getdouble("reservoir_fraction");
    pressure_source = ini_file.getdouble("pressure_source");
    pressure_drain = ini_file.getdouble("pressure_drain");

    dt = ini_file.getdouble("dt");
}
