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
    collide();
}

void System::move() {
    timer->start_moving();
    for(int n=0;n<thread_control.num_molecules;n++) {
        mover->move_molecule(n,dt,rnd,0);
    }
    thread_control.update_molecule_cells_local();
    timer->end_moving();
}

void System::collide() {
    timer->start_colliding();

    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *cell = thread_control.cells[i];

        cell->prepare();
        collisions += cell->collide(rnd);
    }

    timer->end_colliding();
}

void System::accelerate() {
    if(settings->gravity_direction < 0) return;
    for(int n=0;n<thread_control.num_molecules;n++) {
        thread_control.v[3*n+settings->gravity_direction] += settings->gravity*dt;
    }
}
