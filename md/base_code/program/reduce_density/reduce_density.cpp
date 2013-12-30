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
	if(args < 3) {
		cout << "Please specify the number of cpus, relative density" << endl;
		cout << "./reduce_density int double" << endl;
		return 0;
	}

	int cpus = atoi(argv[1]);
	double density = atof(argv[2]);

	double data[6*MAX_ATOM_NUM];
	unsigned long *atom_type = new unsigned long[MAX_ATOM_NUM];
	char *filename = new char[100];
	int total_num_atoms = 0;
	int num_particles;
	for(int cpu=0;cpu<cpus;cpu++) { 
		// Read binary file
		sprintf(filename,"state_files/state%04d.bin",cpu);
		ifstream state_file(filename,ios::in | ios::binary);
		state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
		state_file.read(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
		state_file.read(reinterpret_cast<char*>(atom_type),num_particles*sizeof(unsigned long));
		total_num_atoms += num_particles;
	}

	int num_atoms_to_be_removed = total_num_atoms*(1 - density);
	int num_atoms_to_be_removed0 = num_atoms_to_be_removed;
	while(num_atoms_to_be_removed > 0) {
		for(int cpu=0;cpu<cpus;cpu++) { 
			if(num_atoms_to_be_removed == 0) continue;

			sprintf(filename,"state_files/state%04d.bin",cpu);
			ifstream state_file(filename,ios::in | ios::binary);
			state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
			state_file.read(reinterpret_cast<char*>(&data),6*num_particles*sizeof(double));
			state_file.read(reinterpret_cast<char*>(atom_type),num_particles*sizeof(unsigned long));

			// Loop through all atoms and mark them after the criteria
			int num_conserved_particles = 0;
			for(int i=0;i<num_particles;i++) {
				if(num_atoms_to_be_removed == 0 || rnd(0,1) < density) {
					// Keep this particle
					data[6*num_conserved_particles + 0] = data[6*i + 0];
					data[6*num_conserved_particles + 1] = data[6*i + 1];
					data[6*num_conserved_particles + 2] = data[6*i + 2];
					data[6*num_conserved_particles + 3] = data[6*i + 3];
					data[6*num_conserved_particles + 4] = data[6*i + 4];
					data[6*num_conserved_particles + 5] = data[6*i + 5];
					atom_type[num_conserved_particles] = atom_type[num_conserved_particles];

					num_conserved_particles++;
				} else num_atoms_to_be_removed--; // Count that we skipped this particle.
			}
			state_file.close();

			ofstream save_state_file(filename,ios::out | ios::binary);
			save_state_file.write(reinterpret_cast<char*>(&num_conserved_particles),sizeof(int));
			save_state_file.write(reinterpret_cast<char*>(&data),6*num_conserved_particles*sizeof(double));
			save_state_file.write(reinterpret_cast<char*>(atom_type),num_conserved_particles*sizeof(unsigned long));
			save_state_file.close();
			cout << "CPU " << cpu << " had " << num_particles << ". I now have " << num_conserved_particles << endl;
		}
	}

	//printf("Sphere created at (%f,%f,%f) with radius %f\n",x,y,z,radius);
	cout << endl << "Removed " << num_atoms_to_be_removed0 << " atoms. The density is now " << density << endl;

	return 0;
}