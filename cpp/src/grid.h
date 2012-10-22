#pragma once
#include <CVector.h>
#include <armadillo>
#include <vector>
#include <system.h>

using namespace arma;

class GridPoint
{
public:
    bool is_wall;
    CVector normal;
    GridPoint(bool wall) { is_wall = wall; normal = CVector(); }
};

class Grid
{
public:
    int cols;
    int rows;
    System *system;
    vector<GridPoint> points;

    Grid(mat M);
    GridPoint *get_grid_point(const int i, const int j);
    GridPoint *get_grid_point(const double x, const double y);
    void calculate_normals();
};
