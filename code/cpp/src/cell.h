#pragma once

class System;
#include <random.h>
#include "molecule.h"
#include "system.h"
#include <armadillo>
#include <CVector.h>

using namespace arma;
using namespace std;

/*
struct col_pairs_greater
{
    bool operator()( const Cell& c1, const Cell& c2 ) const {
            return c1.collision_pairs < c2.collision_pairs;
        }
};
*/

class Cell {
public:
	double volume;
    int pixels; // Used to calculate volume in a cell
    int total_pixels;

	double vr_max;
    System *system;
    double density;
    int average_over;

    int collision_pairs;
    int i,j;
    int particles;
    int particle_capacity;
    unsigned int *particle_indices;


    Cell(System *system);
	void reset();
    int prepare();
    void resize(int n);
    int collide(Random *rnd);

    void update_volume();

    static bool cmp(Cell *c1, Cell *c2);

    /*
    bool operator < (const Cell &c1, const Cell &c2) const {
        return c1.collision_pairs < c2.collision_pairs;
    }

    bool operator > (const Cell &c1, const Cell &c2) const {
        return c1.collision_pairs > c2.collision_pairs;
    }
    */
};
