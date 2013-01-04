#include <sorter.h>
#include <algorithm>
#include <vector>

Sorter::Sorter(System *_system) {
    system = _system;
}

void Sorter::sort_system() {
	//* Find the cell address for each particle

    int N = system->N;
    double width = system->width;
    double height = system->height;
    int cells_x = system->cells_x;
    int cells_y = system->cells_y;

    int i, j;

	//* Count the number of particles in each cell
    for(i=0;i<cells_x;i++)
        for(j=0;j<cells_y;j++)
            system->cells[i][j]->reset();

    int *cell_x = new int[N];
    int *cell_y = new int[N];

    // Loop through and calculate cell index for each particle
	for(int n=0; n<N; n++ ) {
        Molecule *m = system->molecules[n];
        if(!m->active) continue;

        i = (int)((m->r(0)/width)*cells_x);
        j = (int)((m->r(1)/height)*cells_y);

        i = min(max(i,0),cells_x-1); // ensure we are within the limits
        j = min(max(j,0),cells_y-1);

        cell_x[n] = i;
        cell_y[n] = j;

        system->cells[i][j]->particles++;
	}

    // Resize the particle index list
    Cell *c;
    for(i=0;i<cells_x;i++)
        for(j=0;j<cells_y;j++) {
            c = system->cells[i][j];
            if(c->particle_capacity < c->particles)
                c->resize(2*c->particles);
            c->particles = 0; // Will use to track which particles we have added
        }

    // Update particle index list
    for(int n=0; n<N; n++ ) {
        Molecule *m = system->molecules[n];
        if(!m->active) continue;

        i = cell_x[n];
        j = cell_y[n];
        c = system->cells[i][j];
        c->particle_indices[c->particles++] = n;
    }
    int collisions = 0;

    for(int i=0;i<cells_x;i++)
        for(int j=0;j<cells_y;j++)
            collisions += system->cells[i][j]->prepare();

    delete [] cell_x;
    delete [] cell_y;

    int collisions_per_thread = collisions/4;

    int load_balance_index = 0;
    for(int n=0;n<system->threads;n++)
        system->cells_in_list[n] = 0;

    int number_of_distributed_collisions = 0;
    for(int i=0;i<system->cells_x;i++)
        for(int j=0;j<system->cells_y;j++) {
            // Put this cell in the list for the correct thread
            system->load_balanced_cell_list[load_balance_index][system->cells_in_list[load_balance_index]++] = system->cells[i][j];
            number_of_distributed_collisions += system->cells[i][j]->collision_pairs;

            if(number_of_distributed_collisions > collisions_per_thread*(1+load_balance_index) && load_balance_index < system->threads-1)
                load_balance_index++;
        }
}
