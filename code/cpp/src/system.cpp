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
    // if(myid==0) cout << steps << endl;

    steps += 1;
    t += dt;
    accelerate();
    move();
    collide();
    int num_p = 0;
    MPI_Allreduce(&thread_control.num_particles,&num_p,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    int num_local = 0;
    for(int i=0;i<thread_control.cells.size();i++) {
        num_local += thread_control.cells[i]->molecules.size();
    }

    // cout << myid << " has " << num_local << " particles. Globally: " << num_p << ". And tmp: " << thread_control.tmp_molecules.size() <<  endl;

    // if(myid==0) cout << steps << endl;
}

void System::move() {
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
    thread_control.update_mpi();
}

void System::collide() {
    for(int i=0;i<thread_control.cells.size();i++) {
        thread_control.cells[i]->prepare();
        collisions += thread_control.cells[i]->collide(rnd);
    }
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
