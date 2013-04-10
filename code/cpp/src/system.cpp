#include <system.h>

#include <math.h>
#include <fstream>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <unitconverter.h>
#include <dsmc_io.h>

#include <time.h>
#include <system.inc.cpp>
#include <dsmctimer.h>

void System::step() {
    steps += 1;
    t += dt;
    accelerate();
    move();
    collide();
    maintain_pressure();
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

    for(int i=0;i<thread_control.my_cells.size();i++) {
        Cell *cell = thread_control.my_cells[i];

        cell->prepare();
        collisions += cell->collide(rnd);
    }

    timer->end_colliding();
}

void System::accelerate() {
    if(settings->gravity_direction < 0) return;
    timer->start_accelerate();

    int gravity_dir = settings->gravity_direction;
    double gravity = settings->gravity*dt;

    for(int n=0;n<thread_control.num_molecules;n++) {
        thread_control.v[3*n+gravity_dir] += gravity;
    }
    timer->end_accelerate();
}

void System::find_position_in_reservoirs(double *r, bool find_position_in_source) {
    bool did_collide = true;
    while(did_collide) {
        r[0] = length[0]*rnd->nextDouble();
        r[1] = length[1]*rnd->nextDouble();
        r[2] = length[2]*rnd->nextDouble();
        if(find_position_in_source)  r[settings->gravity_direction] = reservoir_size*rnd->nextDouble();
        else r[settings->gravity_direction] = (length[settings->gravity_direction] - reservoir_size) + reservoir_size*rnd->nextDouble();

        did_collide = *world_grid->get_voxel(r)>=voxel_type_wall;
    }
}

void System::add_molecule_in_pressure_reservoirs(bool add_in_source) {
    int n = thread_control.num_molecules;

    thread_control.v[3*n+0] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    thread_control.v[3*n+1] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    thread_control.v[3*n+2] = rnd->nextGauss()*sqrt(temperature/settings->mass);

    find_position_in_reservoirs(&thread_control.r[3*n],add_in_source);
    thread_control.r0[3*n+0] = thread_control.r[3*n+0];
    thread_control.r0[3*n+1] = thread_control.r[3*n+1];
    thread_control.r0[3*n+2] = thread_control.r[3*n+2];
    thread_control.molecule_moved[n] = false;
    Cell *cell = thread_control.all_cells[thread_control.cell_index_from_position(&thread_control.r[3*n])];
    cell->add_molecule(n,thread_control.molecule_index_in_cell,thread_control.molecule_cell_index);
    thread_control.num_molecules++;
}

bool System::remove_molecule_in_pressure_reservoir(bool remove_from_source) {
    Cell *cell = NULL;

    if(remove_from_source) cell = thread_control.source_reservoir_cells[ thread_control.source_reservoir_cells.size()*rnd->nextDouble() ];
    else cell = thread_control.drain_reservoir_cells[ thread_control.drain_reservoir_cells.size()*rnd->nextDouble() ];

    if(cell->num_molecules>0) {
        // Remove this random molecule
        int molecule_index_in_cell = cell->num_molecules*rnd->nextDouble();
        int molecule_index = cell->molecules[molecule_index_in_cell];
        cell->remove_molecule(molecule_index,thread_control.molecule_index_in_cell);

        // Move the last molecule into that memory location
        int last_molecule_index = thread_control.num_molecules-1;
        int last_molecule_cell_index = thread_control.molecule_cell_index[last_molecule_index];
        int last_molecule_index_in_cell = thread_control.molecule_index_in_cell[last_molecule_index];
        memcpy(&thread_control.r[3*molecule_index],&thread_control.r[3*last_molecule_index],3*sizeof(double));
        memcpy(&thread_control.v[3*molecule_index],&thread_control.v[3*last_molecule_index],3*sizeof(double));
        memcpy(&thread_control.r0[3*molecule_index],&thread_control.r0[3*last_molecule_index],3*sizeof(double));

        cell = thread_control.all_cells[last_molecule_cell_index];
        thread_control.molecule_cell_index[molecule_index] = last_molecule_cell_index;
        thread_control.molecule_index_in_cell[molecule_index] = last_molecule_index_in_cell;
        cell->molecules[last_molecule_index_in_cell] = molecule_index;
        thread_control.num_molecules--;
        return true;
    } else return false;
}

void System::maintain_pressure() {
    long num_molecules_in_source = 0;
    double volume_in_source = 0;
    double pressure_in_source = 0;
    for(int i=0;i<thread_control.source_reservoir_cells.size();i++) {
        Cell *cell = thread_control.source_reservoir_cells[i];
        num_molecules_in_source += cell->num_molecules;
        volume_in_source += cell->volume;
    }
    double temp_over_volume = 0;
    if(volume_in_source>0) {
        int added_molecules = 0;
        temp_over_volume = temperature/volume_in_source;
        pressure_in_source = eff_num*num_molecules_in_source*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_source);
        if(!(steps%100)) cout << "Pressure is=" << pressure_in_source << " wanted pressure: " << wanted_pressure << endl;

        if(pressure_in_source<wanted_pressure) {
            while(pressure_in_source<wanted_pressure) {
                add_molecule_in_pressure_reservoirs(true);
                num_molecules_in_source++;
                pressure_in_source = eff_num*num_molecules_in_source*temp_over_volume;
                added_molecules++;
            }
        } else {
            while(pressure_in_source>settings->pressure_source) {
                if(remove_molecule_in_pressure_reservoir(true)) {
                    num_molecules_in_source--;
                    pressure_in_source = num_molecules_in_source*temp_over_volume;
                    added_molecules--;
                }
            }
        }
        if(!(steps%100)) cout << "Added " << added_molecules << " molecules." << endl;
    }

}
