#include <iostream>
#include <math.h>
#include <fstream>
#include <molecule.h>
#include <system.h>
#include <time.h>
#include <defines.h>
#include "system.inc.cpp"

void System::step() {
    time_t t0;

    steps += 1;
    t += dt;

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

}

void System::move() {
    for(int n=0; n<N; n++ ) {
        molecules[n]->move(dt,rnd);
    }
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell

    // Cell **local_cells = load_balanced_cell_list[thread_id];

    for(int i=0;i<settings->cells_x;i++) {
        for(int j=0;j<settings->cells_y;j++) {
            col += cells[i][j]->collide(rnd);
        }
    }

    // for(int i=0;i<cells_in_list[thread_id];i++)
    //    local_col += local_cells[i]->collide(rnd);

	return col;
}

void System::accelerate() {
    for(int n=0;n<N;n++) {
        if(molecules[n]->r(0) < max_x_acceleration) {
            molecules[n]->v(0) += acceleration*dt;
        }
    }
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
