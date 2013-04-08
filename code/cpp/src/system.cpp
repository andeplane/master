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
#include <dsmctimer.h>

void System::step() {
    steps += 1;
    t += dt;
    accelerate();
    move();
    timer->start_colliding();
    collide();
    timer->end_colliding();
}

void System::move() {
    timer->start_moving();
    for(int n=0;n<thread_control.num_molecules;n++) {
        mover->move_molecule(n,dt,rnd,0);
    }
    double *r = thread_control.r;

    for(int n=0;n<thread_control.num_molecules;n++) {
        int cell_index_new = thread_control.cell_index_from_position(&r[3*n]);
        int cell_index_old = thread_control.molecule_cell_index[n];

        DummyCell *dummy_cell_new = thread_control.dummy_cells[cell_index_new];
        DummyCell *dummy_cell_old = thread_control.dummy_cells[cell_index_old];

        if(cell_index_new != cell_index_old) {
            // We changed cell
            if(thread_control.dummy_cells[cell_index_new]->node_id != myid) {
                // If this is another node, send it to that list
//                int node_id = thread_control.dummy_cells[cidx]->node_id;
//                thread_control.nodes_new_atoms_list[node_id].push_back(molecule);
            } else {
                // If it is this node, just add it to the dummy cell list
                dummy_cell_old->real_cell->remove_molecule(n,thread_control.molecule_index_in_cell);
                dummy_cell_new->real_cell->add_molecule(n,thread_control.molecule_index_in_cell,thread_control.molecule_cell_index);
            }
        }
    }
    timer->end_moving();

//    timer->start_mpi();
//    thread_control.update_mpi();
//    timer->end_mpi();
}

void System::collide() {
    for(int i=0;i<thread_control.cells.size();i++) {
        thread_control.cells[i]->prepare();
        Cell *cell = thread_control.cells[i];
        collisions += cell->collide(rnd);
    }
}

void System::accelerate() {
    if(settings->gravity_direction < 0) return;
    for(int n=0;n<thread_control.num_molecules;n++) {
        thread_control.v[3*n+settings->gravity_direction] += settings->gravity*dt;
    }
}
