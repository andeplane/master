#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <system.h>
#include <settings.h>
#include <mpi.h>
#include <mdio.h>
#include <unitconverter.h>
#include <mdtimer.h>
#include <random.h>
#include <atom.h>
#include <statisticssampler.h>
#include <atom_types.h>

#include <potential_lennard_jones.h>

using namespace std;
System::System() {

}

void System::setup(int myid_, Settings *settings_) {
    mdtimer = new MDTimer();
    mdtimer->start_system_initialize();
    unit_converter = new UnitConverter();

    myid = myid_;
    settings = settings_;
    num_atoms_local = 0;
    num_atoms_free_global = 0;
    num_atoms_all_global = 0;
    num_atoms_ghost = 0;
    num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;

    steps = 0;
    rnd = new Random(-(myid+1));

    init_parameters();

    mdio = new MDIO();
    mdio->setup(this);
    set_topology();
    if(settings->load_state) mdio->load_state_from_file_binary();
    else create_FCC();

    rearrange_frozen_atoms();

    MPI_Allreduce(&num_atoms_local,&num_atoms_all_global,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&num_atoms_free,&num_atoms_free_global,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    mpi_copy();
    calculate_accelerations();
    half_kick();

    if(myid==0) cout << "System size: " << unit_converter->length_to_SI(system_length[0])*1e10 << " Å " << unit_converter->length_to_SI(system_length[1])*1e10 << " Å " << unit_converter->length_to_SI(system_length[2])*1e10 << " Å" << endl;
    if(myid==0) cout << "Atoms: " << num_atoms_all_global << endl;
    if(myid==0) cout << "Free atoms: " << num_atoms_free_global << endl;
    if(myid==0) cout << "Processors: " << num_nodes << " (" << settings->nodes_x << "," << settings->nodes_y << "," << settings->nodes_z << ")" << endl;
    if(myid==0) cout << "Free atoms/processor: " << num_atoms_free_global/num_nodes << endl;



    mdtimer->end_system_initialize();
}

void System::rearrange_frozen_atoms() {
    double positions_copy[num_atoms_local][3];
    double velocities_copy[3*num_atoms_local];
    unsigned long atom_type_copy[num_atoms_local];

    for(i=0;i<num_atoms_local;i++) {
        positions_copy[i][0] = positions[i][0];
        positions_copy[i][1] = positions[i][1];
        positions_copy[i][2] = positions[i][2];

        velocities_copy[3*i+0] = velocities[3*i+0];
        velocities_copy[3*i+1] = velocities[3*i+1];
        velocities_copy[3*i+2] = velocities[3*i+2];
        atom_type_copy[i] = atom_type[i];
    }

    num_atoms_frozen = 0;
    num_atoms_free = 0;
    for(i=0;i<num_atoms_local;i++) {
        if(atom_type_copy[i] == FROZEN) {
            int index = (num_atoms_frozen + num_atoms_free);

            velocities[3*index+0] = 0;
            velocities[3*index+1] = 0;
            velocities[3*index+2] = 0;
            positions[index][0] = positions_copy[i][0];
            positions[index][1] = positions_copy[i][1];
            positions[index][2] = positions_copy[i][2];
            atom_type[index] = atom_type_copy[i];
            num_atoms_frozen++;
        }
    }

    for(i=0;i<num_atoms_local;i++) {
        if(atom_type_copy[i] != FROZEN) {
            int index = (num_atoms_frozen + num_atoms_free);
            velocities[3*index + 0] = velocities_copy[3*i+0];
            velocities[3*index + 1] = velocities_copy[3*i+1];
            velocities[3*index + 2] = velocities_copy[3*i+2];
            positions[index][0] = positions_copy[i][0];
            positions[index][1] = positions_copy[i][1];
            positions[index][2] = positions_copy[i][2];
            atom_type[index] = atom_type_copy[i];
            num_atoms_free++;

        }
    }
}

void System::create_FCC() {
    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};

    double r[3];
    double T = unit_converter->temperature_from_SI(settings->temperature);

    bool warning_shown = false;
    for(int x = 0; x < settings->nodes_x*settings->unit_cells_x; x++) {
        for(int y = 0; y < settings->nodes_y*settings->unit_cells_y; y++) {
            for(int z = 0; z < settings->nodes_z*settings->unit_cells_z; z++) {
                for(int k = 0; k < 4; k++) {
                    // Set positions and type

                    r[0] = (x+xCell[k]) * settings->FCC_b/unit_converter->L0_angstrom - origo[0];
                    r[1] = (y+yCell[k]) * settings->FCC_b/unit_converter->L0_angstrom - origo[1];
                    r[2] = (z+zCell[k]) * settings->FCC_b/unit_converter->L0_angstrom - origo[2];
                    bool is_mine = true;
                    for(i=0;i<3;i++) {
                        if(!(r[i] >= 0 && r[i] < node_length[i])) is_mine = false;
                    }

                    if(is_mine) {
                        // if(n++>10) continue;
                        for(i=0;i<3;i++) {
                            positions[num_atoms_local][i] = r[i];
                            initial_positions[3*num_atoms_local+i] = r[i];
                        }
                        velocities[3*num_atoms_local+0] = rnd->nextGauss()*sqrt(T*mass_inverse);
                        velocities[3*num_atoms_local+1] = rnd->nextGauss()*sqrt(T*mass_inverse);
                        velocities[3*num_atoms_local+2] = rnd->nextGauss()*sqrt(T*mass_inverse);
                        atom_type[num_atoms_local] = ARGON;

                        num_atoms_local++;
                        if(!warning_shown && num_atoms_local >= 0.6*MAX_ATOM_NUM) {
                            cout << "                 ### WARNING ###" << endl;
                            cout << "NUMBER OF PARTICLES IS MORE THAN 0.8*MAX_ATOM_NUM" << endl << endl;
                            warning_shown = true;
                        }
                    }
                }
            }
        }
    }
}

void System::init_parameters() {
    mass_inverse = 1.0/settings->mass;
    r_cut = settings->r_cut;
    dt = settings->dt;
    dt_half = dt/2;
    t = 0;

    num_processors[0] = settings->nodes_x;
    num_processors[1] = settings->nodes_y;
    num_processors[2] = settings->nodes_z;

    node_index[0] = myid/(settings->nodes_y*settings->nodes_z);
    node_index[1] = (myid/settings->nodes_z) % settings->nodes_y;
    node_index[2] = myid%settings->nodes_z;

    // Size of this node
    node_length[0] = settings->unit_cells_x*settings->FCC_b/unit_converter->L0_angstrom;
    node_length[1] = settings->unit_cells_y*settings->FCC_b/unit_converter->L0_angstrom;
    node_length[2] = settings->unit_cells_z*settings->FCC_b/unit_converter->L0_angstrom;

    for(a=0;a<3;a++) {
        system_length[a] = node_length[a]*num_processors[a];
        origo[a] = (double)node_index[a] * node_length[a];
        num_cells_local[a] = node_length[a]/r_cut;
        num_cells_including_ghosts[a] = num_cells_local[a]+2;

        cell_length[a] = node_length[a]/num_cells_local[a];
    }
    volume = system_length[0]*system_length[1]*system_length[2];

    num_cells_including_ghosts_yz = num_cells_including_ghosts[1]*num_cells_including_ghosts[2];
    num_cells_including_ghosts_xyz = num_cells_including_ghosts_yz*num_cells_including_ghosts[0];
    head_all_atoms = new int[num_cells_including_ghosts_xyz];
    head_free_atoms = new int[num_cells_including_ghosts_xyz];

    for(int cx=0;cx<num_cells_local[0]+2;cx++) {
        for(int cy=0;cy<num_cells_local[1]+2;cy++) {
            for(int cz=0;cz<num_cells_local[2]+2;cz++) {
                cell_index_from_ijk(cx,cy,cz,cell_index);
                if(cx == 0 || cx == num_cells_local[0]+1 || cy == 0 || cy == num_cells_local[1]+1 || cz == 0 || cz == num_cells_local[2]+1) {
                    is_ghost_cell[cell_index] = true;
                } else is_ghost_cell[cell_index] = false;
            }
        }
    }
}

void System::set_topology() {
    /*----------------------------------------------------------------------
    Defines a logical network topology.  Prepares a neighbor-node ID table,
    nn, & a shift-vector table, sv, for internode message passing.  Also
    prepares the node parity table, myparity.
    ----------------------------------------------------------------------*/

    /* Integer vectors to specify the six neighbor nodes */
    int integer_vector[6][3] = {
        {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}
    };

    int k1[3];

    /* Set up neighbor tables, nn & sv */
    for (int n=0; n<6; n++) {
        /* Vector index of neighbor ku */
        for (a=0; a<3; a++) {
            k1[a] = (node_index[a]+integer_vector[n][a]+num_processors[a])%num_processors[a];
        }

        /* Scalar neighbor ID, nn */
        neighbor_nodes[n] = k1[0]*num_processors[1]*num_processors[2]+k1[1]*num_processors[2]+k1[2];
        /* Shift vector, sv */
        for (a=0; a<3; a++) shift_vector[n][a] = node_length[a]*integer_vector[n][a];
    }



    /* Set up the node parity table, myparity */
    for (a=0; a<3; a++) {
        if (num_processors[a] == 1) {
            my_parity[a] = 2;
        } else if (node_index[a]%2 == 0) {
            my_parity[a] = 0;
        } else {
            my_parity[a] = 1;
        }
    }
}

inline void System::cell_index_from_ijk(const int &i, const int &j, const int &k, unsigned int &cell_index) {
    cell_index = i*num_cells_including_ghosts_yz+j*num_cells_including_ghosts[2]+k;
}

inline void System::cell_index_from_vector(unsigned int *mc, unsigned int &cell_index) {
    cell_index = mc[0]*num_cells_including_ghosts_yz+mc[1]*num_cells_including_ghosts[2]+mc[2];
}

inline bool System::atom_should_be_copied(double* ri, int ku) {
  int dimension,higher;
  dimension = ku/2; /* x(0)|y(1)|z(2) direction */
  higher = ku%2; /* Lower(0)|higher(1) direction */
  if (higher == 0) return ri[dimension] < r_cut;
  else return ri[dimension] > node_length[dimension]-r_cut;
}


inline bool System::atom_did_change_node(double* ri, int ku) {
    int dimension,higher;
    dimension = ku/2;    /* x(0)|y(1)|z(2) direction */
    higher = ku%2; /* Lower(0)|higher(1) direction */
    if (higher == 0) return ri[dimension] < 0.0;
    else return ri[dimension] > node_length[dimension];
}

void System::mpi_move() {
    MPI_Status status;

    int new_atoms = 0;
    int node_id,num_send,num_receive;
    short node_higher, node_lower, local_node_id;

    /* Reset the # of to-be-moved atoms, move_queue[][0] */
    for (short ku=0; ku<6; ku++) move_queue[ku][0] = 0;

    /* Make a moved-atom list, move_queue----------------------------------*/

    for(short dimension=0;dimension<3;dimension++) {
        /* Scan all the residents & immigrants to list moved-out atoms */
        for (i=num_atoms_frozen; i<num_atoms_local+new_atoms; i++) {
            node_lower = 2*dimension;
            node_higher = 2*dimension+1;
            /* Register a to-be-copied atom in move_queue[kul|kuh][] */
            if (!atom_moved[i]) { /* Don't scan moved-out atoms */
                // Check if this atom moved
                if (atom_did_change_node(positions[i],node_lower)) move_queue[node_lower][ ++move_queue[node_lower][0] ] = i;
                else if (atom_did_change_node(positions[i],node_higher)) move_queue[node_higher][ ++move_queue[node_higher][0] ] = i;
            }
        }

        /* Message passing with neighbor nodes----------------------------*/

        mdtimer->start_mpi();

        /* Loop over the lower & higher directions------------------------*/
        for (j=0; j<2; j++) {
            local_node_id=2*dimension+j;
            node_id = neighbor_nodes[local_node_id]; /* Neighbor node ID */
            /* Send atom-number information---------------------------------*/

            num_send = move_queue[local_node_id][0]; /* # of atoms to-be-sent */

            /* Even node: send & recv*/
            /* Odd node: recv & send */
            if (my_parity[dimension] == 0) {
                MPI_Send(&num_send,1,MPI_INT,node_id,110,MPI_COMM_WORLD);
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,110,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,110, MPI_COMM_WORLD,&status);
                MPI_Send(&num_send,1,MPI_INT,node_id,110,MPI_COMM_WORLD);
            }
            else num_receive = num_send;

            /* Message buffering */
            for (i=1; i<=num_send; i++) {
                /* Shift the coordinate origin */
                mpi_send_buffer[6*(i-1)    + 0] = positions[ move_queue[local_node_id][i]][0] - shift_vector[local_node_id][0];
                mpi_send_buffer[6*(i-1)    + 1] = positions[ move_queue[local_node_id][i]][1] - shift_vector[local_node_id][1];
                mpi_send_buffer[6*(i-1)    + 2] = positions[ move_queue[local_node_id][i]][2] - shift_vector[local_node_id][2];
                mpi_send_buffer[6*(i-1)+ 3 + 0] = velocities[3*move_queue[local_node_id][i]+0];
                mpi_send_buffer[6*(i-1)+ 3 + 1] = velocities[3*move_queue[local_node_id][i]+1];
                mpi_send_buffer[6*(i-1)+ 3 + 2] = velocities[3*move_queue[local_node_id][i]+2];
                atom_moved[ move_queue[local_node_id][i] ] = true;
            }

            /* Even node: send & recv, if not empty */
            /* Odd node: recv & send, if not empty */
            if (my_parity[dimension] == 0) {
                MPI_Send(mpi_send_buffer,6*num_send,MPI_DOUBLE,node_id,120,MPI_COMM_WORLD);
                MPI_Recv(mpi_receive_buffer,6*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,120,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(mpi_receive_buffer,6*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,120,MPI_COMM_WORLD,&status);
                MPI_Send(mpi_send_buffer,6*num_send,MPI_DOUBLE,node_id,120,MPI_COMM_WORLD);
            }
            /* Single layer: Exchange information with myself */
            else memcpy(mpi_receive_buffer,mpi_send_buffer,6*num_receive*sizeof(double));
            //else for (i=0; i<6*num_receive; i++) mpi_receive_buffer[i] = mpi_send_buffer[i];


            /* Message storing */
            for (i=0; i<num_receive; i++) {
                atom_type[num_atoms_local+new_atoms+i] = ARGON;
                atom_moved[num_atoms_local+new_atoms+i] = false;
                positions [(num_atoms_local+new_atoms+i)][0] = mpi_receive_buffer[6*i   + 0];
                positions [(num_atoms_local+new_atoms+i)][1] = mpi_receive_buffer[6*i   + 1];
                positions [(num_atoms_local+new_atoms+i)][2] = mpi_receive_buffer[6*i   + 2];

                velocities[3*(num_atoms_local+new_atoms+i) + 0] = mpi_receive_buffer[6*i+3 + 0];
                velocities[3*(num_atoms_local+new_atoms+i) + 1] = mpi_receive_buffer[6*i+3 + 1];
                velocities[3*(num_atoms_local+new_atoms+i) + 2] = mpi_receive_buffer[6*i+3 + 2];
            }

            /* Increment the # of new atoms */
            new_atoms += num_receive;

            /* Internode synchronization */
            MPI_Barrier(MPI_COMM_WORLD);
        } /* Endfor lower & higher directions, j */

        mdtimer->end_mpi();
    }

    int ipt = num_atoms_frozen;
    for (i=num_atoms_frozen; i<num_atoms_local+new_atoms; i++) {
        if (!atom_moved[i]) {
            for (a=0; a<3; a++) {
                positions [ipt][a] = positions [i][a];
                velocities[3*ipt+a] = velocities[3*i+a];
            }
            atom_type[ipt] = atom_type[i];

            ipt++;
        }
    }

    /* Update the compressed # of resident atoms */
    num_atoms_local = ipt;
    num_atoms_free = num_atoms_local - num_atoms_free;
}

void System::mpi_copy() {
    MPI_Status status;
    int node_id, num_send, num_receive;
    int new_ghost_atoms = 0;
    short higher, local_node_id;
    for(short dimension=0;dimension<3;dimension++) {
        for (higher=0; higher<2; higher++) move_queue[2*dimension+higher][0] = 0;
        for(i=0;i<num_atoms_local+new_ghost_atoms;i++) {
            for(higher=0;higher<2;higher++) {
                local_node_id = 2*dimension + higher;
                if (atom_should_be_copied(positions[i],local_node_id)) move_queue[local_node_id][++(move_queue[local_node_id][0])] = i;
            }
        }

        mdtimer->start_mpi();

        /* Loop through higher and lower node in this dimension */
        for(higher=0;higher<2;higher++) {
            local_node_id= 2*dimension+higher;
            node_id = neighbor_nodes[local_node_id];
            num_send = move_queue[local_node_id][0];

            if (my_parity[dimension] == 0) {
                MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
                MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);
            }
            else num_receive = num_send;

            for (i=1; i<=num_send; i++) {
                /* Shift the coordinate origin */
                mpi_send_buffer[3*(i-1)+0] = positions[ move_queue[local_node_id][i]][0]-shift_vector[local_node_id][0];
                mpi_send_buffer[3*(i-1)+1] = positions[ move_queue[local_node_id][i]][1]-shift_vector[local_node_id][1];
                mpi_send_buffer[3*(i-1)+2] = positions[ move_queue[local_node_id][i]][2]-shift_vector[local_node_id][2];
            }

            if (my_parity[dimension] == 0) {
                MPI_Send(mpi_send_buffer,3*num_send,MPI_DOUBLE,node_id,20,MPI_COMM_WORLD);
                MPI_Recv(mpi_receive_buffer,3*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
            }
            else if (my_parity[dimension] == 1) {
                MPI_Recv(mpi_receive_buffer,3*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
                MPI_Send(mpi_send_buffer,3*num_send,MPI_DOUBLE,node_id,20,MPI_COMM_WORLD);
            }
            else memcpy(mpi_receive_buffer,mpi_send_buffer,3*num_receive*sizeof(double));
            // else for (i=0; i<3*num_receive; i++) mpi_receive_buffer[i] = mpi_send_buffer[i];

            for (i=0; i<num_receive; i++) {
                positions[ (num_atoms_local+new_ghost_atoms+i) ][0] = mpi_receive_buffer[3*i+0];
                positions[ (num_atoms_local+new_ghost_atoms+i) ][1] = mpi_receive_buffer[3*i+1];
                positions[ (num_atoms_local+new_ghost_atoms+i) ][2] = mpi_receive_buffer[3*i+2];
            }

            new_ghost_atoms += num_receive;
            MPI_Barrier(MPI_COMM_WORLD);
        }

        mdtimer->end_mpi();
    }

    num_atoms_ghost = new_ghost_atoms;
}

void System::half_kick() {
    for(n=num_atoms_frozen;n<num_atoms_local;n++) {
        velocities[3*n+0] += accelerations[3*n+0]*dt_half;
        velocities[3*n+1] += accelerations[3*n+1]*dt_half;
        velocities[3*n+2] += accelerations[3*n+2]*dt_half;
    }
}

void System::full_kick() {
    for(n=num_atoms_frozen;n<num_atoms_local;n++) {
        velocities[3*n+0] += accelerations[3*n+0]*dt;
        velocities[3*n+1] += accelerations[3*n+1]*dt;
        velocities[3*n+2] += accelerations[3*n+2]*dt;
    }
}

void System::move() {
    mdtimer->start_moving();

    for(n=num_atoms_frozen;n<num_atoms_local;n++) {
        positions[n][0] += velocities[3*n+0]*dt;
        positions[n][1] += velocities[3*n+1]*dt;
        positions[n][2] += velocities[3*n+2]*dt;
        atom_moved[n] = false;
    }

    mdtimer->end_moving();
}

void System::apply_gravity() {
    for(n=num_atoms_frozen;n<num_atoms_local;n++) {
        velocities[3*n+settings->gravity_direction] += settings->gravity_force*dt;
    }
}

void System::step() {
    move();
    mpi_move();
    mpi_copy();
    if(settings->gravity_direction >= 0) apply_gravity();
    if(settings->many_frozen_atoms) calculate_accelerations_many_frozen_atoms();
    else calculate_accelerations();
    full_kick();
    steps++;
    t += dt;
}
