#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define MAX_ATOM_NUM 100000
#define ARGON  0
#define FROZEN 1
double pos_coeff = 3.405;
double pos_coeff_squared = pos_coeff*pos_coeff;

using namespace std;

int main(int args, char *argv[]) {
	if(args < 7) {
		cout << "Please specify the number of cpus, radius, cylinder dimension, center_x, center_y, center_z" << endl;
		cout << "./create_cylinder int double [0|1|2] double double double" << endl;
		return 0;
	}

	int cpus = atoi(argv[1]);
	double radius = atof(argv[2]);
	int dimension = atoi(argv[3]);
	double center[3];
	center[0] = atof(argv[4])/pos_coeff;
	center[1] = atof(argv[5])/pos_coeff;
	center[2] = atof(argv[6])/pos_coeff;

	double data[6*MAX_ATOM_NUM];
	unsigned long  atom_type[MAX_ATOM_NUM];
	char *filename = new char[100];
	int num_particles;
	
	double radius_squared = radius*radius;
	double r2;
	int frozen_atoms = 0;
	for(int cpu=0;cpu<cpus;cpu++) { 
		sprintf(filename,"state_files/state%04d.bin",cpu);
		ifstream state_file(filename,ios::in | ios::binary);

		state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
		state_file.read(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
		state_file.read(reinterpret_cast<char*>(&atom_type),num_particles*sizeof(unsigned long));
		
		for(int i=0;i<num_particles;i++) {
			r2 = 0;
			for(int a=0;a<3;a++) {
				if(a != dimension) r2 += (data[6*i+a]-center[a])*(data[6*i+a]-center[a])*pos_coeff_squared;
			}
			
			if(r2 > radius_squared) {
				atom_type[i] = FROZEN;
				frozen_atoms++;
			}
		}
		state_file.close();

		ofstream save_state_file(filename,ios::out | ios::binary);
		save_state_file.write(reinterpret_cast<char*>(&num_particles),sizeof(int));
		save_state_file.write(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
		save_state_file.write(reinterpret_cast<char*>(&atom_type),num_particles*sizeof(unsigned long));
		save_state_file.close();
		cout << "Created cylinder in cpu " << cpu << endl;
	}

	cout << "Cylinder created with " << frozen_atoms << " frozen atoms." << endl;
	
	return 0;
}