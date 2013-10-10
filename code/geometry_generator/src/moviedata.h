#pragma once
#include <vector>
#include <string>
#include <defines.h>
using std::vector;
using std::string;

class Timestep
{
public:
    Timestep();
    Timestep *next;
    Timestep *previous;
    vector<double> positions;
    vector<double> velocities;
    void add_molecule_data(vector<float> &positions);
#ifdef OPENGL
    void render();
#endif
};

class MovieData
{
public:
    int cpus;
    int timesteps;
    Timestep *first_timestep;
    MovieData(int cpus_, int timesteps_);
    void load_movie_files(string state_folder);
};
