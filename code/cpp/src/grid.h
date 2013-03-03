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
    bool is_wall_boundary;
    int i,j;

    double T;
    CVector normal;
    CVector tangent;
    GridPoint(bool wall) { is_wall = wall;
                           is_wall_boundary = false;
                           normal = CVector();
                           tangent = CVector();
                           T = 3.0; }
};

class Grid
{
public:
    int cols;
    int rows;
    System *system;
    vector<GridPoint> points;

    Grid(mat &M, System *_system);
    inline GridPoint *get_grid_point(const int &i, const int &j);
    GridPoint *get_grid_point(const double &x, const double &y);
    GridPoint* get_grid_point(double *r);
    void calculate_inner_points();
    void calculate_normals();
    void calculate_tangents();
};
