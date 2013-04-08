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
    thread_control.update_molecule_cells_local();
    timer->end_moving();

//    for(int dimension = 0; dimension<3; dimension++) {
//        timer->start_moving();
//        thread_control.update_molecule_cells(dimension);
//        thread_control.update_new_molecules(dimension);
//        timer->end_moving();

//        timer->start_mpi();
//        thread_control.update_mpi(dimension);
//        timer->end_mpi();
//    }

//    thread_control.update_molecule_arrays();
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
