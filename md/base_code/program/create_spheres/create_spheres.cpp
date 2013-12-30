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

inline double rnd(double min, double max) {
	return min + (double)rand()/RAND_MAX*(max-min);
}

int main(int args, char *argv[]) {
	if(args < 8) {
		cout << "Please specify the number of cpus, num_spheres, r_min, r_max, x_max, y_max, z_max" << endl;
		cout << "./create_cylinder int int double double   double double double   " << endl;
		return 0;
	}

	int cpus = atoi(argv[1]);
	int num_spheres = atoi(argv[2]);
	double r_min = atof(argv[3])/pos_coeff;
	double r_max = atof(argv[4])/pos_coeff;
	double x_max = atof(argv[5])/pos_coeff;
	double y_max = atof(argv[6])/pos_coeff;
	double z_max = atof(argv[7])/pos_coeff;
	
	double x_max_half = x_max/2;
	double y_max_half = y_max/2;
	double z_max_half = z_max/2;

	double data[6*MAX_ATOM_NUM];
	unsigned long  *atom_type = new unsigned long[MAX_ATOM_NUM];
	char *filename = new char[100];
	int frozen_atoms = 0;
	double max_x = 0;
	
	for(int n=0;n<num_spheres;n++) {
		double radius = rnd(r_min,r_max);
		double radius_squared = radius*radius;
		double x = rnd(0,x_max);
		double y = rnd(0,y_max);
		double z = rnd(0,z_max);
		
		int num_particles;
		for(int cpu=0;cpu<cpus;cpu++) { 
			// Read binary file
			sprintf(filename,"state_files/state%04d.bin",cpu);
			ifstream state_file(filename,ios::in | ios::binary);
			state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
			state_file.read(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
			state_file.read(reinterpret_cast<char*>(atom_type),num_particles*sizeof(unsigned long));
			
			// Loop through all atoms and mark them after the criteria
			for(int i=0;i<num_particles;i++) {
				double dx = data[6*i + 0] - x;
				double dy = data[6*i + 1] - y;
				double dz = data[6*i + 2] - z;

				// Add periodic boundary conditions
				if(dx > x_max_half) dx -= x_max;
				if(dy > y_max_half) dy -= x_max;
				if(dz > z_max_half) dz -= x_max;
				if(dx < -x_max_half) dx += x_max;
				if(dy < -y_max_half) dy += x_max;
				if(dz < -z_max_half) dz += x_max;

				double dr2 = dx*dx + dy*dy + dz*dz;
				
				if(dr2 < radius_squared) {
					atom_type[i] = FROZEN;
					frozen_atoms++;
				}
			}
			state_file.close();

			ofstream save_state_file(filename,ios::out | ios::binary);
			save_state_file.write(reinterpret_cast<char*>(&num_particles),sizeof(int));
			save_state_file.write(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
			save_state_file.write(reinterpret_cast<char*>(atom_type),num_particles*sizeof(unsigned long));
			save_state_file.close();
		}

		printf("Sphere created at (%f,%f,%f) with radius %f\n",x*pos_coeff,y*pos_coeff,z*pos_coeff,radius*pos_coeff);
		//printf("Sphere created at (%f,%f,%f) with radius %f\n",x,y,z,radius);
	}
	cout << endl << "We now have " << frozen_atoms << " frozen atoms." << endl;

	return 0;
}