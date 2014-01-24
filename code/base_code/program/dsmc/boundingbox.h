#pragma once
#include <vector>
#include <system.h>
using std::vector;

class BoundingBox
{
public:
    System *system;
    double corner[3];
    double size[3];
    vector<int> particle_indices;

    void add_particle(int particle_index);
    BoundingBox(int particle_index, System *system_);
    bool point_is_inside_box(double *r);
    bool did_collide_with_a_particle(double *r);
};
