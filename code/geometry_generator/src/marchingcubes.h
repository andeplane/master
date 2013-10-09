#pragma once
#include <mesh.h>
#include <iostream>
#include <cvector.h>
using std::cout;
using std::endl;

class ComplexGeometry;

typedef struct {
   CVector p[3];
} TRIANGLE;

typedef struct {
   CVector p[8];
   double val[8];
} GRIDCELL;

class MarchingCubes : public Mesh
{
private:
    CVector interpolate_vertices(double isolevel, CVector p1, CVector p2, double valp1, double valp2);
    int polygonise(GRIDCELL grid, double isolevel, TRIANGLE *triangles);

public:
    int num_triangles;
    MarchingCubes():Mesh() { num_triangles = 0;}
    void create_marching_cubes_from_complex_geometry(ComplexGeometry &cg, CVector system_length, const int min_value);
    void load_from_file(string filename);
    void save_to_file(string filename);
};
