#pragma once
#include <settings.h>
#include <system.h>

class DSMC_IO
{
public:
    Settings *settings;
    System *system;

    DSMC_IO(System *system_);
    void save_state_to_file_binary();
    void save_state_to_file_xyz();
    void load_state_from_file_binary();
    void save_state_to_movie_file();
    void finalize();
    bool movie_file_open;
    int  movie_frames;

    FILE *energy_file;
    FILE *movie_file;
};
