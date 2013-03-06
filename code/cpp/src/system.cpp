#include <system.h>

#include <math.h>
#include <fstream>
#include "molecule.h"
#include "cell.h"
#include "sorter.h"
#include <grid.h>
#include "random.h"
#include <unitconverter.h>
#include <dsmc_io.h>

#include <time.h>
#include <defines.h>
#include <system.inc.cpp>


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

    for(int i=0;i<settings->cells_x;i++) {
        for(int j=0;j<settings->cells_y;j++) {
            for(int k=0;k<settings->cells_z;k++) {
                col += cells[i][j][k]->collide(rnd);
            }
        }
    }

	return col;
}

void System::accelerate() {
    for(int n=0;n<N;n++) {
        if(molecules[n]->r[0] < max_x_acceleration) {
            molecules[n]->v[0] += acceleration*dt;
        }
    }
}
