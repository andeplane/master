/*
#include <sorter.h>

#include <algorithm>
#include <vector>
#include <cell.h>
#include <system.h>
#include <molecule.h>

Sorter::Sorter(System *_system) {
    system = _system;
}

void Sorter::sort_system() {
	//* Find the cell address for each particle

    int N = system->N;
    double Lx = system->Lx;
    double Ly = system->Ly;
    double Lz = system->Lz;
    int cells_x = system->settings->cells_x;
    int cells_y = system->settings->cells_y;
    int cells_z = system->settings->cells_z;

    int i, j, k;

	//* Count the number of particles in each cell
    for(i=0;i<cells_x;i++) {
        for(j=0;j<cells_y;j++) {
            for(k=0;k<cells_z;k++) {
                system->cells[i][j][k]->reset();
            }
        }
    }

    int *cell_x = new int[N];
    int *cell_y = new int[N];
    int *cell_z = new int[N];

    // Loop through and calculate cell index for each particle
	for(int n=0; n<N; n++ ) {
        Molecule *m = system->molecules[n];
        if(!m->active) continue;

        i = (int)((m->r[0]/Lx)*cells_x);
        j = (int)((m->r[1]/Ly)*cells_y);
        k = (int)((m->r[2]/Lz)*cells_z);

        i = min(max(i,0),cells_x-1); // ensure we are within the limits
        j = min(max(j,0),cells_y-1);
        k = min(max(k,0),cells_z-1);

        cell_x[n] = i;
        cell_y[n] = j;
        cell_z[n] = k;

        system->cells[i][j][k]->particles++;
	}

    // Resize the particle index list
    Cell *c;
    for(i=0;i<cells_x;i++) {
        for(j=0;j<cells_y;j++) {
            for(k=0;k<cells_z;k++) {
                c = system->cells[i][j][k];
                if(c->particle_capacity < c->particles) {
                    c->resize(2*c->particles);
                }
                c->particles = 0; // Will use to track which particles we have added
            }
        }
    }

    // Update particle index list
    for(int n=0; n<N; n++ ) {
        Molecule *m = system->molecules[n];
        if(!m->active) continue;

        i = cell_x[n];
        j = cell_y[n];
        k = cell_z[n];
        c = system->cells[i][j][k];
        c->particle_indices[c->particles++] = n;
    }
    int collisions = 0;

    for(int i=0;i<cells_x;i++) {
        for(int j=0;j<cells_y;j++) {
            for(int k=0;k<cells_z;k++) {
                collisions += system->cells[i][j][k]->prepare();
            }
        }
    }

    delete [] cell_x;
    delete [] cell_y;
    delete [] cell_z;
}
*/
