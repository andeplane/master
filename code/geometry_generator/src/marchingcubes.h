#pragma once
#include <mesh.h>
#include <iostream>
#include <cvector.h>
using std::cout;
using std::endl;

// #define OLDWAY

class ComplexGeometry;
#ifdef OLDWAY
typedef struct {
   CVector p[3];
} TRIANGLE;

typedef struct {
   CVector p[8];
   double val[8];
} GRIDCELL;
#endif
class MarchingCubes : public Mesh
{
private:
#ifdef OLDWAY
    CVector interpolate_vertices(double isolevel, CVector p1, CVector p2, double valp1, double valp2);
    int polygonise(GRIDCELL grid, double isolevel, TRIANGLE *triangles);
#endif

public:
    int num_triangles;
    MarchingCubes():Mesh() { num_triangles = 0;}
    void create_marching_cubes_from_complex_geometry(ComplexGeometry &cg, CVector system_length, double threshold);
    void load_from_file(string filename);
    void save_to_file(string filename);
};
