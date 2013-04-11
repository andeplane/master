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
    if(myid==0) accelerate();
    move();
    if(myid==0) collide();
    if(myid==0 && settings->maintain_pressure) maintain_pressure();
}

void System::send_molecules_to_slaves() {

    int molecules_per_node = num_molecules/num_nodes;
    for(int node_id=1;node_id<num_nodes;node_id++) {
        int start_index = molecules_per_node*node_id;
        int end_index = start_index + molecules_per_node;
        if(node_id == num_nodes-1) end_index = num_molecules - 1;
        end_index = min(end_index,num_molecules - 1);

        int num_send = end_index - start_index + 1;

        MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);
        MPI_Send(&r[3*start_index],3*num_send,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD);
        MPI_Send(&v[3*start_index],3*num_send,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD);
        MPI_Send(&r0[3*start_index],3*num_send,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD);
    }

    num_molecules_this_node = molecules_per_node;
}

void System::receive_molecules_from_master() {
    MPI_Status status;

    MPI_Recv(&num_molecules,1,MPI_INT,0,10,MPI_COMM_WORLD,&status);
    MPI_Recv(r,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status);
    MPI_Recv(v,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status);
    MPI_Recv(r0,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status);
    num_molecules_this_node = num_molecules;
}

void System::send_molecules_to_master() {
    MPI_Send(r,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
    MPI_Send(v,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
    MPI_Send(r0,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
}

void System::receive_molecules_from_slaves() {
    MPI_Status status;

    int molecules_per_node = num_molecules/num_nodes;
    for(int node_id=1;node_id<num_nodes;node_id++) {
        int start_index = molecules_per_node*node_id;
        int end_index = start_index + molecules_per_node;
        if(node_id == num_nodes-1) end_index = num_molecules - 1;
        end_index = min(end_index,num_molecules - 1);

        int num_receive = end_index - start_index + 1;

        MPI_Recv(&r[3*start_index],3*num_receive,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD,&status);
        MPI_Recv(&v[3*start_index],3*num_receive,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD,&status);
        MPI_Recv(&r0[3*start_index],3*num_receive,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD,&status);
    }
}

void System::move() {
    timer->start_mpi();
    if(myid==0) send_molecules_to_slaves();
    else receive_molecules_from_master();
    timer->end_mpi();

    timer->start_moving();
    // cout << myid << " will move " << num_molecules_this_node << " molecules." << endl;
    for(int n=0;n<num_molecules_this_node;n++) {
        mover->move_molecule(n,dt,rnd,0);
    }
    timer->end_moving();

    timer->start_mpi();
    if(myid==0) receive_molecules_from_slaves();
    else send_molecules_to_master();
    timer->end_mpi();

    timer->start_moving();
    if(myid==0) update_molecule_cells();
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
    double kinetic_energy = 0;

    for(int i=0;i<source_reservoir_cells.size();i++) {
        Cell *cell = source_reservoir_cells[i];
        num_molecules_in_source += cell->num_molecules;
        volume_in_source += cell->volume;
        kinetic_energy += cell->calculate_kinetic_energy();
    }
    double temperature_in_reservoir = 2.0/3*kinetic_energy/num_molecules_in_source;

    double temp_over_volume = 0;
    if(volume_in_source>0) {
        int added_molecules = 0;
        temp_over_volume = temperature_in_reservoir/volume_in_source;
        pressure_in_source = atoms_per_molecule*num_molecules_in_source*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_source);
        long wanted_num_molecules = wanted_pressure*volume_in_source/temperature/atoms_per_molecule;
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
    double kinetic_energy = 0;

    for(int i=0;i<drain_reservoir_cells.size();i++) {
        Cell *cell = drain_reservoir_cells[i];
        num_molecules_in_drain += cell->num_molecules;
        volume_in_drain += cell->volume;
        kinetic_energy += cell->calculate_kinetic_energy();
    }
    double temperature_in_reservoir = 2.0/3*kinetic_energy/num_molecules_in_drain;

    double temp_over_volume = 0;
    if(volume_in_drain>0) {
        int added_molecules = 0;
        temp_over_volume = temperature_in_reservoir/volume_in_drain;
        pressure_in_drain = atoms_per_molecule*num_molecules_in_drain*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_drain);
        long wanted_num_molecules = wanted_pressure*volume_in_drain/temperature/atoms_per_molecule;
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
