#include <iostream>
#include <math.h>
#include <fstream>
#include <molecule.h>
#include <system.h>
#include <time.h>
#include <omp.h>
#include <defines.h>
#include "system.inc.cpp"

void System::step() {
    time_t t0;

    steps += 1;
    t += 1;

    accelerate();

    //* Move all the particles
    t0 = clock();
    move();
    time_consumption[MOVE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

    //* Sort the particles into cells
    t0 = clock();
    // Arrange all particles in correct cells
    sorter->sort_system();

    time_consumption[SORT] += ((double)clock()-t0)/CLOCKS_PER_SEC;

    t0 = clock();

    collisions += collide();
    time_consumption[COLLIDE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

    t0 = clock();

    double E = 0;
    vec p(3,1);
    double density = 0;

    Cell *c;
    #pragma omp parallel for private(c) num_threads(threads)
    for(int i=0;i<cells_x;i++)
        for(int j=0;j<cells_y;j++) {
            c = cells[i][j];
            c->sampleStatistics();
            // E += c->energy;
            // p += c->momentum;

            // density += c->density;
    }

    time_consumption[SAMPLE] += ((double)clock()-t0)/CLOCKS_PER_SEC;
}

void System::move() {
    #pragma omp parallel num_threads(threads)
    {
        Random *rnd = randoms[omp_get_thread_num()];
        #pragma omp for
        for(int n=0; n<N; n++ ) {
            molecules[n]->move(dt,rnd);
        }
    }
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell
    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        Random *rnd = randoms[thread_id];
        int local_col = 0;
        Cell **local_cells = load_balanced_cell_list[thread_id];

        #pragma omp for
        for(int i=0;i<cells_x;i++)
            for(int j=0;j<cells_y;j++)
                local_col += cells[i][j]->collide(rnd);

        #pragma omp critical
        {
           col += local_col;
        }
    }

	return col;
}

void System::accelerate() {
    for(int n=0;n<N;n++)
        if(molecules[n]->r(0) < max_x_acceleration)
            molecules[n]->v(0) += acceleration*dt;
}

void System::printPositionsToFile(FILE *file) {
    fprintf(file,"%d\n",N);
	fprintf(file,"Random comment that must be here\n");

    for(int n=0;n<N;n++) {
        // fprintf(file,"%s %.10f %.10f 0\n",molecules[n]->type,molecules[n]->r(0),molecules[n]->r(1));
        // We return height - r(1) because system is inverted
        fprintf(file,"%s %.10f %.10f 0\n",molecules[n]->information_carrier ? "H" : "H",molecules[n]->r(0),-molecules[n]->r(1) + height);
    }
	
}
