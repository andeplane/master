#pragma once
#include <vector>
#include <string>
#include <defines.h>
#include <marchingcubes.h>

using std::vector;
using std::string;
class CVector;

class Timestep
{
public:
    Timestep(CVector system_length);
    Timestep *next;
    Timestep *previous;
    vector<float> positions;
    vector<float> velocities;
    CVector system_length;
    bool is_marching_cubes_built;
    void add_molecule_data(vector<float> &positions);
    void scale_positions(float scale_factor);
    MarchingCubes marching_cubes;
#ifdef OPENGL
    void render_points();
    void build_marching_cubes();
    void render_marching_cubes();
#endif
};

class MovieData
{
public:
    int cpus;
    int timesteps;
    Timestep *first_timestep;
    MovieData(int cpus_, int timesteps_);
    void load_movie_files(string state_folder, CVector system_length);
};
