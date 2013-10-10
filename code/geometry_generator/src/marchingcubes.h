#pragma once
#include <mesh.h>
#include <iostream>
#include <cvector.h>
using std::cout;
using std::endl;

class ComplexGeometry;

class MarchingCubes : public Mesh
{
public:
    int num_triangles;
    MarchingCubes():Mesh() { num_triangles = 0;}
    template <typename T>
    void create_marching_cubes_from_array(const T* scalar_field, int nx, int ny, int nz, CVector box_length, double threshold);
    void create_marching_cubes_from_complex_geometry(ComplexGeometry &cg, CVector system_length, double threshold);
    void load_from_file(string filename);
    void save_to_file(string filename);
};
