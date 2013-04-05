#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

void freeze_atoms() {
    int num_frozen = 0;
    for(n=0;n<num_atoms_local;n++) {
        double dr2 = 0;
        double dx;

        for(a=0;a<3;a++) {
            dx = positions[3*n+a] - 7.5;

            dr2 += pow(dx,2);
        }

        if(sqrt(dr2) > 5) {
            cout << "I froze atom " << num_frozen++ << endl;
            frozen_atom[n] = true;
            velocities[3*n+0] = 0;
            velocities[3*n+1] = 0;
            velocities[3*n+2] = 0;
        }

        atom_moved[n] = false;
    }
}

int main(int args, char *argv[]) {
	if(args < 2) {
		cout << "Please specify the number of cpus." << endl;
		return 0;
	}

	int cpus = atoi(argv[1]);
	ifstream **state_files = new ifstream*[cpus];
	double **positions = new double*[1000000*cpus];
	for(int cpu=0;cpu<cpus;cpu++) {
		char *filename = new char[100];
		sprintf(filename,"release/state_files/movie%04d.bin",cpu);
		movie_files[cpu] = new ifstream(filename,ios::in | ios::binary);
	}
	cout << cpus << " state files opened." << endl;

	for(int timestep=0;timestep<timesteps;timestep++) {
		int num_particles = 0;
		for(int cpu=0;cpu<cpus;cpu++) { 
			int N;
			movie_files[cpu]->read(reinterpret_cast<char*>(&N),sizeof(int));
			movie_files[cpu]->read(reinterpret_cast<char*>(&positions[3*num_particles]),3*N*sizeof(double));
			num_particles += N;
		}
		

		file << num_particles << endl;
		file << "sup" << endl;
		for(int n=0;n<num_particles;n++) {
        	// We return height - r(1) because system is inverted
        	file << "H " << positions[3*n+0] << " " << positions[3*n+1] << " " << positions[3*n+2] << endl;
	    }

	    cout << "Wrote timestep " << timestep << endl;
	}

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		movie_files[cpu]->close();
	}
	cout << "Movie created." << endl;

	return 0;
}