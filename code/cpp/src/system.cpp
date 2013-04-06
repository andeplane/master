#include <system.h>

#include <math.h>
#include <fstream>
#include <molecule.h>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <unitconverter.h>
#include <dsmc_io.h>
#include <settings.h>
#include <time.h>
#include <defines.h>
#include <system.inc.cpp>
#include <dsmctimer.h>

void System::step() {
    steps += 1;
    t += dt;
    // accelerate();
    move();
    timer->start_colliding();
    collide();
    timer->end_colliding();
}

void System::move() {
    timer->start_moving();

    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *cell = thread_control.cells[i];
        mover->move_molecules(cell,dt,rnd);
    }

    timer->end_moving();

    for(int dimension = 0;dimension<3;dimension++) {
        for(int i=0;i<thread_control.cells.size();i++) {
            Cell *cell = thread_control.cells[i];
            cell->update_molecule_cells(dimension);
        }
        timer->start_mpi();
        thread_control.update_mpi(dimension);
        timer->end_mpi();
    }

    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *cell = thread_control.cells[i];
        cell->update_molecule_arrays();
    }


}

void System::collide() {
    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *cell = thread_control.cells[i];

        cell->prepare();
        collisions += cell->collide(rnd);
    }
}

void System::accelerate() {
    int k = 0;
    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *c = thread_control.cells[i];
        for(int j=0;j<c->num_molecules;j++) {
            if(c->r[3*j+0] < max_x_acceleration) {
                k++;
                c->v[3*j+0] += acceleration*dt;
            }
        }
    }
}
