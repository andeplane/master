#include <complexgeometry.h>
#include <perlin.h>
#include <progressbar.h>
#include <fstream>
#include <cmath>
#include <random.h>
#include <fstream>
#include <cvector.h>

using std::ifstream;
using std::ofstream;
using std::ios;
using std::fstream;

ComplexGeometry::ComplexGeometry()
{
    vertices = NULL;
    normals = NULL;
    tangents1 = NULL;
    tangents2 = NULL;
    nx = 0; ny = 0; nz = 0; num_vertices = 0;
    has_normals_tangents_and_boundary = false;
    max_value = 0;
}

void ComplexGeometry::allocate(int nx_, int ny_, int nz_) {
    nx = nx_; ny = ny_; nz = nz_;
    num_vertices = nx*ny*nz;

    vertices = new unsigned char[num_vertices];
    vertices_float = new float[num_vertices];
    normals = new float[3*num_vertices];
    tangents1 = new float[3*num_vertices];
    tangents2 = new float[3*num_vertices];
}

void ComplexGeometry::load_text_files(string base_filename, CVector matrix_size, double threshold) {
    matrix_size.z = 2;
    allocate(matrix_size.x, matrix_size.y, matrix_size.z);
    int file_start = 31;
    int num_files = nz;
    char filename[1000];
    int num_vertices_so_far = 0;
    for(int file_num = file_start; file_num<file_start+nz; file_num++) {
        sprintf(filename, "%s%02d.txt",base_filename.c_str(), file_num);
        ifstream infile(filename);
        for(int line_num=0; line_num<nx; line_num++) {
            for(int i=0; i<ny; i++) {
                float value;
                infile >> value;
                max_value = max(max_value,value);
                vertices_float[num_vertices_so_far] = value;
                vertices[num_vertices_so_far++] = value >= threshold;
            }
        }
        infile.close();
    }

    // calculate_normals_tangents_and_inner_points(1.0);
    // save_to_vtk("mus.vtk");
}

void ComplexGeometry::save_to_vtk(string filename) {
    ofstream ofile(filename.c_str());

    int N = nx*ny*nz;

    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << "structured point" << endl;
    ofile << "ASCII" << endl;
    ofile << endl;
    ofile << "DATASET STRUCTURED_POINTS" << endl;
    ofile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    ofile << "ORIGIN 0.0 0.0 0.0" << endl;
    // ofile << "SPACING 1 1 1" << endl;
    ofile << "SPACING " << 1.0/double(nx) << " " << 1.0/double(ny) << " " << 1.0/double(nz) << endl;
    ofile << "POINT_DATA " << N << endl;
    ofile << "SCALARS atomdist double" << endl;
    ofile << "LOOKUP_TABLE default" << endl;
    ofile << endl;

    // column-major ordering...
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int index = i + j*nx + k*nx*ny;
                ofile << vertices[index] << endl;
            }
        }
    }

    ofile.close();
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
    has_normals_tangents_and_boundary = true;
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



void ComplexGeometry::create_perlin_geometry(int nx_, int ny_, int nz_, int octave, int frequency, int amplitude , int seed, float threshold, bool do_calculate_normals_tangents_and_boundary) {
    Perlin p(octave, frequency, amplitude, seed);
    allocate(nx_, ny_, nz_);

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

    if(do_calculate_normals_tangents_and_boundary) {
        calculate_normals_tangents_and_inner_points(1);
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
