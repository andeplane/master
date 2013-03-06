#include <system.h>

#include <math.h>
#include <fstream>
#include <molecule.h>
#include <cell.h>
#include <sorter.h>
#include <grid.h>
#include <random.h>
#include <unitconverter.h>
#include <dsmc_io.h>

#include <time.h>
#include <defines.h>
#include <system.inc.cpp>


void System::step() {
    steps += 1;
    t += dt;
    // accelerate();
    move();
    collisions += collide();
}

void System::move() {
    int cidx;
    Molecule *current_molecule;
    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *c = thread_control.cells[i];
        Molecule *m = c->first_molecule;
        while(m != NULL) {
            current_molecule = m;
            m = m->next;

            current_molecule->move(dt,rnd);
            cidx = thread_control.cell_index_from_molecule(current_molecule);
            if(cidx != current_molecule->cell_index) {
                // We changed cell
                thread_control.dummy_cells[cidx]->new_molecules.push_back(current_molecule);
                thread_control.dummy_cells[m->cell_index]->real_cell->remove_molecule(current_molecule);
            }
        }
    }

    thread_control.update_mpi();
}

int System::collide() {
    /*
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
    */
}

void System::accelerate() {
//    for(int n=0;n<N;n++) {
//        if(molecules[n]->r[0] < max_x_acceleration) {
//            molecules[n]->v[0] += acceleration*dt;
//        }
//    }
}
