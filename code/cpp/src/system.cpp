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
    accelerate();

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
        timer->start_moving();
        for(int i=0;i<thread_control.cells.size();i++) {
            Cell *cell = thread_control.cells[i];
            cell->update_molecule_cells(dimension);
        }
        thread_control.update_new_molecules(dimension);

        timer->end_moving();
        timer->start_mpi();
        thread_control.update_mpi(dimension);
        timer->end_mpi();
    }

    timer->start_moving();
    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *cell = thread_control.cells[i];
        cell->update_molecule_cells_local();
    }

    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *cell = thread_control.cells[i];
        cell->update_molecule_arrays();
    }
    thread_control.add_new_molecules_to_cells();
    timer->end_moving();
}

void System::collide() {
    int molecules = 0;
    int molecules_squared = 0;
    int max = 0;
    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *cell = thread_control.cells[i];
        int mol = cell->num_molecules;
        if(mol>max) max = mol;
        molecules += mol;
        molecules_squared += mol*mol;
        cell->prepare();
        collisions += cell->collide(rnd);
    }
    molecules /= thread_control.cells.size();
    molecules_squared /= thread_control.cells.size();
}

void System::accelerate() {
    if(settings->gravity_direction < 0) return;

    for(int i=0;i<thread_control.cells.size();i++) {
        Cell *c = thread_control.cells[i];
        for(int j=0;j<c->num_molecules;j++) {
            c->v[3*j+settings->gravity_direction] += settings->gravity*dt;
        }
    }
}
