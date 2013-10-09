#pragma once
#include <string>
#include <vector>
using std::string;
using std::vector;

class ComplexGeometry
{
public:
    unsigned char *vertices;
    float *normals;
    float *tangents1;
    float *tangents2;
    unsigned int nx, ny, nz, num_vertices;
    bool has_normals_tangents_and_boundary;

    ComplexGeometry();

    void load_from_binary_file_without_normals_and_tangents(string filename, bool calculate_normals_and_tangents, int number_of_neighbor_averages);
    void create_perlin_geometry(int nx, int ny, int nz, int octave, int frequency, int amplitude , int seed, float threshold, bool do_calculate_normals_tangents_and_boundary);
    void calculate_normals_tangents_and_inner_points(int number_of_neighbor_averages);
    void find_boundary_points();
    void calculate_tangents();
    void calculate_normals(int number_of_neighbor_average);
    void save_to_file(string filename);
};
