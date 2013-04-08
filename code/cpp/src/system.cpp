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
    int cidx;
    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *c = thread_control.cells[i];
        int size = c->molecules.size();

        for(int j=0;j<size;j++) {
            Molecule *molecule = c->molecules[j];
            molecule->move(dt,rnd);
            cidx = thread_control.cell_index_from_molecule(molecule);
            if(cidx != molecule->cell_index) {
                // We changed cell
                if(thread_control.dummy_cells[cidx]->node_id != myid) {
                    // If this is another node, send it to that list
                    int node_id = thread_control.dummy_cells[cidx]->node_id;
                    thread_control.nodes_new_atoms_list[node_id].push_back(molecule);
                } else {
                    // If it is this node, just add it to the dummy cell list
                    thread_control.dummy_cells[cidx]->new_molecules.push_back(molecule);
                }
            }
        }
    }
    thread_control.update_local_cells();
    timer->end_moving();

    timer->start_mpi();
    thread_control.update_mpi();
    timer->end_mpi();
}

void System::collide() {
    for(int i=0;i<thread_control.cells.size();i++) {
        thread_control.cells[i]->prepare();
        collisions += thread_control.cells[i]->collide(rnd);
    }
}

void System::accelerate() {
    if(settings->gravity_direction < 0) return;

        for(int i=0;i<thread_control.cells.size();i++) {
            Cell *c = thread_control.cells[i];

            for(int j=0;j<c->molecules.size();j++) {
                Molecule *m = c->molecules[j];
                m->v[settings->gravity_direction] += settings->gravity*dt;
            }
        }
}
