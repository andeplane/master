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
    if(settings->maintain_pressure) maintain_pressure();
//    double pressure = num_molecules*eff_num/volume*temperature;
//    cout << "Global pressure: " << unit_converter->pressure_to_SI(pressure) << endl;
}

void System::move() {
    timer->start_moving();
    for(int n=0;n<num_molecules;n++) {
        mover->move_molecule(n,dt,rnd,0);
    }
    update_molecule_cells();
    timer->end_moving();
}

void System::collide() {
    timer->start_colliding();

    for(int i=0;i<active_cells.size();i++) {
        Cell *cell = active_cells[i];

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

    for(int n=0;n<num_molecules;n++) {
        v[3*n+gravity_dir] += gravity;
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
    int n = num_molecules;

    v[3*n+0] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    v[3*n+1] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    v[3*n+2] = rnd->nextGauss()*sqrt(temperature/settings->mass);

    find_position_in_reservoirs(&r[3*n],add_in_source);
    r0[3*n+0] = r[3*n+0];
    r0[3*n+1] = r[3*n+1];
    r0[3*n+2] = r[3*n+2];

    Cell *cell = all_cells[cell_index_from_position(&r[3*n])];
    cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
    num_molecules++;
}

bool System::remove_molecule_in_pressure_reservoir(bool remove_from_source) {
    Cell *cell = NULL;

    if(remove_from_source) cell = source_reservoir_cells[ source_reservoir_cells.size()*rnd->nextDouble() ];
    else cell = drain_reservoir_cells[ drain_reservoir_cells.size()*rnd->nextDouble() ];

    if(cell->num_molecules>0) {
        // Remove this random molecule
        int this_molecule_index_in_cell = cell->num_molecules*rnd->nextDouble();
        int molecule_index = cell->molecules[this_molecule_index_in_cell];
        cell->remove_molecule(molecule_index,molecule_index_in_cell);

        // Move the last molecule into that memory location
        int last_molecule_index = num_molecules-1;

        while(last_molecule_index==molecule_index) {
            last_molecule_index--;
        }

        int last_molecule_cell_index = molecule_cell_index[last_molecule_index];
        int last_molecule_index_in_cell = molecule_index_in_cell[last_molecule_index];
        memcpy(&r[3*molecule_index],&r[3*last_molecule_index],3*sizeof(double));
        memcpy(&v[3*molecule_index],&v[3*last_molecule_index],3*sizeof(double));
        memcpy(&r0[3*molecule_index],&r0[3*last_molecule_index],3*sizeof(double));

        cell = all_cells[last_molecule_cell_index];
        molecule_cell_index[molecule_index] = last_molecule_cell_index;
        molecule_index_in_cell[molecule_index] = last_molecule_index_in_cell;
        cell->molecules[last_molecule_index_in_cell] = molecule_index;
        num_molecules--;
        return true;
    } else return false;
}

void System::update_molecule_cells() {
    for(int n=0;n<num_molecules;n++) {
        int cell_index_new = cell_index_from_position(&r[3*n]);
        int cell_index_old = molecule_cell_index[n];

        Cell *new_cell = all_cells[cell_index_new];
        Cell *old_cell = all_cells[cell_index_old];

        if(cell_index_new != cell_index_old) {
            // We changed cell
            old_cell->remove_molecule(n,molecule_index_in_cell);
            new_cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
        }
    }
}

void System::maintain_pressure() {
    timer->start_pressure();
    maintain_pressure_source();
    maintain_pressure_drain();
    timer->end_pressure();
}

void System::maintain_pressure_source() {
    long num_molecules_in_source = 0;
    double volume_in_source = 0;
    double pressure_in_source = 0;
    for(int i=0;i<source_reservoir_cells.size();i++) {
        Cell *cell = source_reservoir_cells[i];
        num_molecules_in_source += cell->num_molecules;
        volume_in_source += cell->volume;
    }
    double temp_over_volume = 0;
    if(volume_in_source>0) {
        int added_molecules = 0;
        temp_over_volume = temperature/volume_in_source;
        pressure_in_source = eff_num*num_molecules_in_source*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_source);
        long wanted_num_molecules = wanted_pressure*volume_in_source/temperature/eff_num;
        long delta = wanted_num_molecules-num_molecules_in_source;

        if(pressure_in_source<wanted_pressure) {
            for(int i=0;i<abs(delta);i++) {
                add_molecule_in_pressure_reservoirs(true);
                added_molecules++;
            }
        } else {
            for(int i=0;i<abs(delta);i++) {
                if(remove_molecule_in_pressure_reservoir(true)) {
                    added_molecules--;
                } else i--;
            }
        }
    }
}

void System::maintain_pressure_drain() {
    long num_molecules_in_drain = 0;
    double volume_in_drain = 0;
    double pressure_in_drain = 0;
    for(int i=0;i<drain_reservoir_cells.size();i++) {
        Cell *cell = drain_reservoir_cells[i];
        num_molecules_in_drain += cell->num_molecules;
        volume_in_drain += cell->volume;
    }

    double temp_over_volume = 0;
    if(volume_in_drain>0) {
        int added_molecules = 0;
        temp_over_volume = temperature/volume_in_drain;
        pressure_in_drain = eff_num*num_molecules_in_drain*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_drain);
        long wanted_num_molecules = wanted_pressure*volume_in_drain/temperature/eff_num;
        long delta = wanted_num_molecules-num_molecules_in_drain;

        if(pressure_in_drain<wanted_pressure) {
            for(int i=0;i<abs(delta);i++) {
                add_molecule_in_pressure_reservoirs(false);
                added_molecules++;
            }
        } else {
            for(int i=0;i<abs(delta);i++) {
                if(remove_molecule_in_pressure_reservoir(false)) {
                    added_molecules--;
                } else i--;
            }
        }
    }
}
