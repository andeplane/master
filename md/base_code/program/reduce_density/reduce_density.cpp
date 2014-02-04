#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define MAX_ATOM_NUM 10000000
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

	double *data = new double[6*MAX_ATOM_NUM];
	unsigned long *atom_type = new unsigned long[MAX_ATOM_NUM];
	unsigned long *atom_ids = new unsigned long[MAX_ATOM_NUM];
	char *filename = new char[100];
	int total_free_atoms = 0;
	int num_particles;
	for(int cpu=0;cpu<cpus;cpu++) { 
		// Read binary file
		sprintf(filename,"state_files/state%04d.bin",cpu);
		ifstream state_file(filename,ios::in | ios::binary);
		state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
		state_file.read(reinterpret_cast<char*>(data),6*num_particles*sizeof(double));
		state_file.read(reinterpret_cast<char*>(atom_type),num_particles*sizeof(unsigned long));
		state_file.read(reinterpret_cast<char*>(atom_ids),num_particles*sizeof(unsigned long));
		for(int n=0; n<num_particles; n++) {
			if(atom_type[n] != FROZEN) {
				total_free_atoms++;
			}
		}
	}

	int num_atoms_to_be_removed = total_free_atoms*(1 - density);
	int num_atoms_to_be_removed0 = num_atoms_to_be_removed;
	total_free_atoms = 0; // Recount in loop 

	while(num_atoms_to_be_removed > 0) {
		for(int cpu=0;cpu<cpus;cpu++) { 
			if(num_atoms_to_be_removed == 0) continue;

			sprintf(filename,"state_files/state%04d.bin",cpu);
			ifstream state_file(filename,ios::in | ios::binary);
			state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(int));
			state_file.read(reinterpret_cast<char*>(data),6*num_particles*sizeof(double));
			state_file.read(reinterpret_cast<char*>(atom_type),num_particles*sizeof(unsigned long));
			state_file.read(reinterpret_cast<char*>(atom_ids),num_particles*sizeof(unsigned long));

			// Loop through all atoms and mark them after the criteria
			int num_conserved_particles = 0;
			for(int n=0;n<num_particles;n++) {
				if(atom_type[n] == FROZEN || num_atoms_to_be_removed == 0 || rnd(0,1) < density) {
					// Keep this particle
					data[6*num_conserved_particles + 0] = data[6*n + 0];
					data[6*num_conserved_particles + 1] = data[6*n + 1];
					data[6*num_conserved_particles + 2] = data[6*n + 2];
					data[6*num_conserved_particles + 3] = data[6*n + 3];
					data[6*num_conserved_particles + 4] = data[6*n + 4];
					data[6*num_conserved_particles + 5] = data[6*n + 5];
					atom_type[num_conserved_particles] = atom_type[n];
					atom_ids[num_conserved_particles] = atom_ids[n];

					num_conserved_particles++;
					total_free_atoms += (atom_type[n] != FROZEN); // Count new total number of free atoms

				} else num_atoms_to_be_removed--; // Count that we skipped this particle.
			}
			state_file.close();

			ofstream save_state_file(filename,ios::out | ios::binary);
			save_state_file.write(reinterpret_cast<char*>(&num_conserved_particles),sizeof(int));
			save_state_file.write(reinterpret_cast<char*>(data),6*num_conserved_particles*sizeof(double));
			save_state_file.write(reinterpret_cast<char*>(atom_type),num_conserved_particles*sizeof(unsigned long));
			save_state_file.write(reinterpret_cast<char*>(atom_ids),num_conserved_particles*sizeof(unsigned long));
			save_state_file.close();
			cout << "CPU " << cpu << " had " << num_particles << ". I now have " << num_conserved_particles << endl;
		}
	}

	//printf("Sphere created at (%f,%f,%f) with radius %f\n",x,y,z,radius);

	cout << endl << "Removed " << num_atoms_to_be_removed0 << " atoms. The density is now " << density << endl;
	ofstream free_atoms_file("number_of_free_atoms.txt");
	free_atoms_file << total_free_atoms;
	free_atoms_file.close();

	return 0;
}