#include <complexgeometry.h>
#include <perlin.h>
#include <progressbar.h>
#include <fstream>
#include <cmath>
#include <random.h>

using std::ifstream;
using std::ofstream;
using std::ios;

ComplexGeometry::ComplexGeometry()
{

}

void ComplexGeometry::load_from_binary_file_without_normals_and_tangents(string filename, bool calculate_normals_and_tangents, int number_of_neighbor_averages) {
    ifstream file (filename.c_str(), ios::in | ios::binary);
    file.read (reinterpret_cast<char*>(&nx), sizeof(unsigned int));
    file.read (reinterpret_cast<char*>(&ny), sizeof(unsigned int));
    file.read (reinterpret_cast<char*>(&nz), sizeof(unsigned int));
    num_vertices = nx*ny*nz;
    vertices = new unsigned char[num_vertices];

    file.read (reinterpret_cast<char*>(vertices), num_vertices*sizeof(unsigned char));
    file.close();

    if(calculate_normals_and_tangents) calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::calculate_normals_tangents_and_inner_points(int number_of_neighbor_averages) {
    calculate_normals(number_of_neighbor_averages);
    calculate_tangents();
    find_boundary_points();
}

void ComplexGeometry::save_to_file(string filename) {
    ofstream file (filename.c_str(), ios::out | ios::binary);
    file.write (reinterpret_cast<char*>(&nx), sizeof(unsigned int));
    file.write (reinterpret_cast<char*>(&ny), sizeof(unsigned int));
    file.write (reinterpret_cast<char*>(&nz), sizeof(unsigned int));
    file.write (reinterpret_cast<char*>(vertices), num_vertices*sizeof(unsigned char));
    file.write (reinterpret_cast<char*>(normals),   3*num_vertices*sizeof(float));
    file.write (reinterpret_cast<char*>(tangents1), 3*num_vertices*sizeof(float));
    file.write (reinterpret_cast<char*>(tangents2), 3*num_vertices*sizeof(float));
    file.close();
}



void ComplexGeometry::create_perlin_geometry(int nx_, int ny_, int nz_, int octave, int frequency, int amplitude , int seed, float threshold) {
    Perlin p(octave, frequency, amplitude, seed);
    nx = nx_; ny = ny_; nz = nz_; num_vertices = nx*ny*nz;

    vertices = new unsigned char[num_vertices];
    normals = new float[3*num_vertices];
    tangents1 = new float[3*num_vertices];
    tangents2 = new float[3*num_vertices];

    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                double x = (i-nx/2.0)/(double)nx;
                double y = (j-ny/2.0)/(double)ny;
                double z = (k-nz/2.0)/(double)nz;
                float s = 1.0;

                int index = i + j*nx + k*nx*ny;
                double val = 0;
                for (int a=0; a<5  ; a++) {
                    s = 3.0*a + 2.2513531;

                    val += p.Get(x*s, y*s, z*s);
                }
                // val = pow(val,4.0)*cos(val);;
                if(val < threshold) vertices[index] = 0;
                else vertices[index] = 1;
            }
        }
    }
}

void ComplexGeometry::calculate_normals(int number_of_neighbor_average) {
    bool at_least_one_wall_neighbor;
    bool all_neighbors_are_walls;
    ProgressBar progress_bar(nx, "Creating normals");
    for(int i=0;i<nx;i++) {
        progress_bar.update(i);
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                at_least_one_wall_neighbor = false;
                all_neighbors_are_walls = true;

                int idx = i + j*nx + k*nx*ny;
                for(int di=-1;di<=1;di++) {
                    for(int dj=-1;dj<=1;dj++) {
                        for(int dk=-1;dk<=1;dk++) {
                            int idx2 = ((i+di+nx)%nx) + ((j+dj+ny)%ny)*nx+ ((k+dk+nz)%nz)*nx*ny;

                            if(vertices[idx2]>0) {
                                // If at least one wall neighbor, this is not a single wall voxel
                                at_least_one_wall_neighbor = true;
                            }
                            if(vertices[idx2]==0) {
                                // If not all neighbors are walls and we have norm=1, this is a
                                // single plane that has no defined normal vector.
                                all_neighbors_are_walls = false;
                            }

                            normals[3*idx+0] -= vertices[idx2]*di;
                            normals[3*idx+1] -= vertices[idx2]*dj;
                            normals[3*idx+2] -= vertices[idx2]*dk;
                        }
                    }
                }

                double norm_squared = normals[3*idx+0]*normals[3*idx+0] + normals[3*idx+1]*normals[3*idx+1] + normals[3*idx+2]*normals[3*idx+2];

                if(norm_squared > 0) {
                    normals[3*idx+0] = 0;
                    normals[3*idx+1] = 0;
                    normals[3*idx+2] = 0;
                    for(int di=-number_of_neighbor_average; di<=number_of_neighbor_average; di++) {
                        for(int dj=-number_of_neighbor_average; dj<=number_of_neighbor_average; dj++) {
                            for(int dk=-number_of_neighbor_average; dk<=number_of_neighbor_average; dk++) {
                                int idx2 = ((i+di+nx)%nx) + ((j+dj+ny)%ny)*nx+ ((k+dk+nz)%nz)*nx*ny;

                                normals[3*idx+0] -= vertices[idx2]*di;
                                normals[3*idx+1] -= vertices[idx2]*dj;
                                normals[3*idx+2] -= vertices[idx2]*dk;
                            }
                        }
                    }
                }


                double norm = sqrt(normals[3*idx+0]*normals[3*idx+0] + normals[3*idx+1]*normals[3*idx+1] + normals[3*idx+2]*normals[3*idx+2]);

                if(norm > 0) {
                    normals[3*idx+0] /= norm;
                    normals[3*idx+1] /= norm;
                    normals[3*idx+2] /= norm;
                } else if(!at_least_one_wall_neighbor || !all_neighbors_are_walls) {
                    // Single point or single pixel-plane, should not be a wall
                    vertices[idx] = 0;
                }
            }
        }
    }
}

void ComplexGeometry::calculate_tangents() {
    Random *rnd = new Random(-1);
    ProgressBar progress_bar(nx, "Creating tangent vectors");
    for(int i=0;i<nx;i++) {
        progress_bar.update(i);
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int idx = i + j*nx + k*nx*ny;
                tangents1[3*idx+0] = rnd->next_double();
                tangents1[3*idx+1] = rnd->next_double();
                tangents1[3*idx+2] = rnd->next_double();

                double dot_product = normals[3*idx+0]*tangents1[3*idx+0] + normals[3*idx+1]*tangents1[3*idx+1] + normals[3*idx+2]*tangents1[3*idx+2];

                // Perform gram-schmidt
                tangents1[3*idx+0] -= normals[3*idx+0]*dot_product;
                tangents1[3*idx+1] -= normals[3*idx+1]*dot_product;
                tangents1[3*idx+2] -= normals[3*idx+2]*dot_product;

                // Normalize
                double norm = sqrt(tangents1[3*idx+0]*tangents1[3*idx+0] + tangents1[3*idx+1]*tangents1[3*idx+1] + tangents1[3*idx+2]*tangents1[3*idx+2]);

                if(norm>0) {
                    tangents1[3*idx+0] /= norm;
                    tangents1[3*idx+1] /= norm;
                    tangents1[3*idx+2] /= norm;
                }

                // t2 = n x t1
                tangents2[3*idx+0] = tangents1[3*idx+1]*normals[3*idx+2] - tangents1[3*idx+2]*normals[3*idx+1];
                tangents2[3*idx+1] = tangents1[3*idx+2]*normals[3*idx+0] - tangents1[3*idx+0]*normals[3*idx+2];
                tangents2[3*idx+2] = tangents1[3*idx+0]*normals[3*idx+1] - tangents1[3*idx+1]*normals[3*idx+0];

                // Normalize
                norm = sqrt(tangents2[3*idx+0]*tangents2[3*idx+0] + tangents2[3*idx+1]*tangents2[3*idx+1] + tangents2[3*idx+2]*tangents2[3*idx+2]);

                if(norm>0) {
                    tangents2[3*idx+0] /= norm;
                    tangents2[3*idx+1] /= norm;
                    tangents2[3*idx+2] /= norm;
                }
            }
        }
    }
}

void ComplexGeometry::find_boundary_points() {
    ProgressBar progress_bar(nx, "Finding boundary points");
    for(int i=0;i<nx;i++) {
        progress_bar.update(i);
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int idx = i + j*nx + k*nx*ny;
                double normal_norm = normals[3*idx+0]*normals[3*idx+0] + normals[3*idx+1]*normals[3*idx+1] + normals[3*idx+2]*normals[3*idx+2];

                if(vertices[idx] > 0 && normal_norm>0) {
                    vertices[idx] = 2;
                }
            }
        }
    }
}
