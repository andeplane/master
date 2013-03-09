#pragma once
#include <stdio.h>
#include <fstream>

using namespace std;

class Settings;
class System;

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
    double *positions;

    FILE *energy_file;
    ofstream *movie_file;
};
