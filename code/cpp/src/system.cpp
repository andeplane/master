#include <system.h>

#include <math.h>
#include <fstream>
#include <molecule.h>
#include <cell.h>
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
    accelerate();
    move();
    collisions += collide();
    cout << steps << endl;
}

void System::move() {
    int cidx;
    double sum_pos = 0;
    double sum_pos_after = 0;

    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *c = thread_control.cells[i];
        int size = c->molecules.size();

        for(int j=0;j<size;j++) {
            Molecule *molecule = c->molecules[j];
            molecule->move(dt,rnd);

            // sum_pos += norm(molecule->r,2);

            cidx = thread_control.cell_index_from_molecule(molecule);
            if(cidx != molecule->cell_index) {
                // We changed cell
                if(thread_control.dummy_cells[cidx]->node_id != myid) {
                    // If this is another node, send it to that list
                    thread_control.nodes_new_atoms_list[thread_control.dummy_cells[cidx]->node_id].push_back(molecule);
                } else {
                    // If it is this node, just add it to the dummy cell list
                    thread_control.dummy_cells[cidx]->new_molecules.push_back(molecule);
                }
            }
        }
    }
    thread_control.update_local_cells();
    thread_control.update_mpi();

//    for(int i=0;i<thread_control.cells.size();i++) {
//        Cell *c = thread_control.cells[i];
//        for(int j=0;j<c->molecules.size();j++) {
//            Molecule *molecule = c->molecules[j];
//            sum_pos_after += norm(molecule->r,2);
//        }
//    }

    // cout << "Sum pos before: " << sum_pos << endl;
    // cout << "Sum pos after: " << sum_pos_after << endl;
}

int System::collide() {
	int col = 0;          // Count number of collisions

    for(int i=0;i<thread_control.cells.size();i++) {
        col += thread_control.cells[i]->collide(rnd);
    }

    return col;
}

void System::accelerate() {
    // if(steps>100) return;
    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *c = thread_control.cells[i];
        for(int j=0;j<c->molecules.size();j++) {
            Molecule *molecule = c->molecules[j];
            if(molecule->r[0] < max_x_acceleration) {
                molecule->v[0] += acceleration*dt;
            }
        }
    }
}
