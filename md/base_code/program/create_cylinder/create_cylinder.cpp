#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define MAX_ATOM_NUM 1000000
#define ARGON  0
#define FROZEN 1
double pos_coeff = 3.405;

using namespace std;

int main(int args, char *argv[]) {
	if(args < 7) {
		cout << "Please specify the number of cpus, radius, number of cylinders per dimension, lx, ly" << endl;
		cout << "./create_cylinder int double int double double" << endl;
		return 0;
	}

	int cpus = atoi(argv[1]);
	double radius = atof(argv[2]);
	int num_cylinders_per_dimension = atoi(argv[3]);
	double system_length[2];
	system_length[0] = atof(argv[4]);
	system_length[1] = atof(argv[5]);
	int num_cylinders = num_cylinders_per_dimension*num_cylinders_per_dimension;

	cout << "Will create " << num_cylinders << " cylinders with radii " << radius << " on a system of size (" << system_length[0] << "," << system_length[1] << ")" << endl;
	
	double data[6*MAX_ATOM_NUM];
	unsigned long  atom_type[MAX_ATOM_NUM];
	char *filename = new char[100];
	int num_particles;
	
	double radius_squared = radius*radius;
	int frozen_atoms = 0;

	double cylinder_center_displacement_x = system_length[0] / num_cylinders_per_dimension;
	double cylinder_center_displacement_y = system_length[1] / num_cylinders_per_dimension;
	
	for(int cpu=0;cpu<cpus;cpu++) { 
		sprintf(filename,"state_files/state%04d.bin",cpu);
		ifstream state_file(filename,ios::in | ios::binary);

		state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
		state_file.read(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
		state_file.read(reinterpret_cast<char*>(&atom_type),num_particles*sizeof(unsigned long));
		
		for(int i=0;i<num_particles;i++) {
			double x = data[6*i+0]*pos_coeff;
			double y = data[6*i+1]*pos_coeff;
			bool freed_this_atom = false;
			
			atom_type[i] = FROZEN;
			frozen_atoms++;

			for(int cylinder_x=0; cylinder_x<num_cylinders_per_dimension; cylinder_x++) {
                for(int cylinder_y=0; cylinder_y<num_cylinders_per_dimension; cylinder_y++) {
                    double cylinder_center_x = cylinder_center_displacement_x*(cylinder_x + 0.5);
                    double cylinder_center_y = cylinder_center_displacement_y*(cylinder_y + 0.5);
                    double dx = (x - cylinder_center_x);
                    double dy = (y - cylinder_center_y);
                    double dr2 = dx*dx + dy*dy;

                    if(dr2 < radius_squared) {
                        atom_type[i] = ARGON;
                        frozen_atoms--;
						freed_this_atom = true;
                    }

                    if(freed_this_atom) break;
                }

                if(freed_this_atom) break;
            }
		}
		state_file.close();

		ofstream save_state_file(filename,ios::out | ios::binary);
		save_state_file.write(reinterpret_cast<char*>(&num_particles),sizeof(int));
		save_state_file.write(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
		save_state_file.write(reinterpret_cast<char*>(&atom_type),num_particles*sizeof(unsigned long));
		save_state_file.close();
	}
	cout << num_cylinders << " cylinders created with " << frozen_atoms << " frozen atoms." << endl;
	
	return 0;
}