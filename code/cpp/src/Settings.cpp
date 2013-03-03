#include "settings.h"

Settings::Settings(string filename) {
    ini_file.load(filename);

    load_previous_state = ini_file.getbool("load_previous_state");
    create_movie = ini_file.getbool("create_movie");
    movie_every_n_frame = ini_file.getint("movie_every_n_frame");
    number_of_particles = ini_file.getint("N");
    timesteps = ini_file.getint("timesteps");
    temperature = ini_file.getdouble("temperature");
    wall_temperature = ini_file.getdouble("wall_temperature");
    dt_factor = ini_file.getdouble("dt_factor");
    acceleration = ini_file.getdouble("acceleration");
    max_x_acceleration = ini_file.getdouble("max_x_acceleration");
    density = ini_file.getdouble("density");
    diam = ini_file.getdouble("diam");
    width = ini_file.getdouble("width");
    height = ini_file.getdouble("height");
    cells_x = ini_file.getint("cells_x");
    cells_y = ini_file.getint("cells_y");
}
