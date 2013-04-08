#include <threadcontrol.h>
#include <system.h>
#include <settings.h>
#include <cell.h>
#include <molecule.h>
#include <random.h>
#include <grid.h>
#include <mpi.h>
#include <dsmc_io.h>
#include <structs.h>

inline int sign(int a) { return (a>0) ? 1 : -1; }

void ThreadControl::setup(System *system_) {
    system = system_;
    settings = system->settings;
    num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;
    myid = system->myid;

    my_vector_index[0] = myid/(settings->nodes_y*settings->nodes_z); // Node id in x-direction
    my_vector_index[1] = (myid/settings->nodes_z) % settings->nodes_y; // Node id in y-direction
    my_vector_index[2] = myid % settings->nodes_z; // Node id in z-direction
    num_processors[0] = settings->nodes_x;
    num_processors[1] = settings->nodes_y;
    num_processors[2] = settings->nodes_z;

    setup_topology();

    for(int a=0;a<3;a++) {
        origo[a] = my_vector_index[a]*(system->length[a]/num_processors[a]);
    }

    if(myid==0) cout << "Setting up cells..." << endl;
    setup_cells();
    if(myid==0) cout << "Setting up molecules..." << endl;
    setup_molecules();
    mpi_receive_buffer= new double[9*MAX_PARTICLE_NUM];

    molecules_to_be_moved = new double*[6];
    num_molecules_to_be_moved = new int[6];
    for(int i=0;i<6;i++) {
        molecules_to_be_moved[i] = new double[9*MAX_PARTICLE_NUM];
        num_molecules_to_be_moved[i] = 0;
    }

    new_molecules = new double[9*MAX_PARTICLE_NUM];
    molecule_moved = new bool[MAX_PARTICLE_NUM];
    num_new_molecules = 0;
}

void ThreadControl::setup_topology() {
    int integer_vector[6][3] = {
            {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}
        };

    int k1[3];

    /* Set up neighbor tables */
    for (int n=0; n<6; n++) {
        /* Vector index of neighbor */
        for (int a=0; a<3; a++) {
            k1[a] = (my_vector_index[a]+integer_vector[n][a]+num_processors[a])%num_processors[a];
        }

        /* Scalar neighbor ID, nn */
        neighbor_nodes[n] = k1[0]*num_processors[1]*num_processors[2]+k1[1]*num_processors[2]+k1[2];
    }
}

void ThreadControl::update_cell_volume() {
    for(int i=0;i<cells.size();i++) {
        Cell *c = cells[i];
        c->vr_max = 3*system->mpv;
        c->update_volume();
    }
}

void ThreadControl::setup_molecules() {
    num_particles = settings->number_of_particles*porosity/num_nodes;
    if(settings->load_previous_state) {
        system->io->load_state_from_file_binary();
        return;
    }

    double r[3];
    double v[3];

    for(int n=0;n<num_particles;n++) {
        find_position(r);
        v[0] = system->rnd->nextGauss()*sqrt(system->temperature);
        v[1] = system->rnd->nextGauss()*sqrt(system->temperature);
        v[2] = system->rnd->nextGauss()*sqrt(system->temperature);

        dummy_cells[cell_index_from_position(r)]->real_cell->add_molecule(r,v);
    }
}

inline void ThreadControl::find_position(double *r) {
    bool did_collide = true;
    while(did_collide) {
        r[0] = origo[0] + system->Lx/settings->nodes_x*system->rnd->nextDouble();
        r[1] = origo[1] + system->Ly/settings->nodes_y*system->rnd->nextDouble();
        r[2] = origo[2] + system->Lz/settings->nodes_z*system->rnd->nextDouble();

        did_collide = *system->world_grid->get_voxel(r)>=voxel_type_wall;
    }
}

inline int ThreadControl::cell_index_from_ijk(const int &i, const int &j, const int &k) {
    return i*system->cells_y*system->cells_z + j*system->cells_z + k;
}

int ThreadControl::cell_index_from_position(double *r) {
    int i = r[0]/system->Lx*system->cells_x;
    int j = r[1]/system->Ly*system->cells_y;
    int k = r[2]/system->Lz*system->cells_z;

    return cell_index_from_ijk(i,j,k);
}

int relative_node_index(int index0, int index1, int cells) {
    if(index0 != index1 && index0 == 0 && index1 == cells-1) return -1;
    if(index0 != index1 && index0 == cells-1 && index1 == 0) return 1;
    else return index1-index0;
}

void ThreadControl::setup_cells() {
    int num_cells_global = num_nodes*settings->cells_per_node_x*settings->cells_per_node_y*settings->cells_per_node_z;
    dummy_cells.reserve(num_cells_global);

    for(int i=0;i<system->cells_x;i++) {
        for(int j=0;j<system->cells_y;j++) {
            for(int k=0;k<system->cells_z;k++) {

                int cell_node_id_x = i/settings->cells_per_node_x;
                int cell_node_id_y = j/settings->cells_per_node_y;
                int cell_node_id_z = k/settings->cells_per_node_z;
                int node_id = cell_node_id_x*settings->nodes_y*settings->nodes_z + cell_node_id_y*settings->nodes_z + cell_node_id_z;

                DummyCell *dc = new DummyCell();
                dc->node_index_vector[0] = cell_node_id_x; dc->node_index_vector[1] = cell_node_id_y; dc->node_index_vector[2] = cell_node_id_z;
                dc->node_delta_index_vector[0] = relative_node_index(my_vector_index[0],cell_node_id_x,system->cells_x);
                dc->node_delta_index_vector[1] = relative_node_index(my_vector_index[1],cell_node_id_y,system->cells_y);
                dc->node_delta_index_vector[2] = relative_node_index(my_vector_index[2],cell_node_id_z,system->cells_z);


                dc->node_id = node_id;
                dc->index = cell_index_from_ijk(i,j,k);
                dummy_cells.push_back(dc);

                if(node_id == myid) {
                    Cell *c = new Cell(system);
                    c->index = cell_index_from_ijk(i,j,k);
                    c->vr_max = 3*system->mpv;
                    cells.push_back(c);
                    dc->real_cell = c;
                    c->dummy_cell = dc;
                }
            }
        }
    }
    calculate_porosity();
}

void ThreadControl::calculate_porosity() {
    int filled_pixels = 0;
    int all_pixels = 0;

    int i_start = float(my_vector_index[0])*settings->cells_per_node_x/system->cells_x*system->world_grid->Nx;
    int i_end   = float(my_vector_index[0]+1)*settings->cells_per_node_x/system->cells_x*system->world_grid->Nx;

    int j_start = float(my_vector_index[1])*settings->cells_per_node_y/system->cells_y*system->world_grid->Ny;
    int j_end   = float(my_vector_index[1]+1)*settings->cells_per_node_y/system->cells_y*system->world_grid->Ny;

    int k_start = float(my_vector_index[2])*settings->cells_per_node_z/system->cells_z*system->world_grid->Nz;
    int k_end   = float(my_vector_index[2]+1)*settings->cells_per_node_z/system->cells_z*system->world_grid->Nz;
    int cell_index, c_x, c_y, c_z;

    for(int k=k_start;k<k_end;k++) {
        for(int j=j_start;j<j_end;j++) {
            for(int i=i_start;i<i_end;i++) {
                c_x = i*system->cells_x/(float)system->world_grid->Nx;
                c_y = j*system->cells_y/(float)system->world_grid->Ny;
                c_z = k*system->cells_z/(float)system->world_grid->Nz;
                cell_index = cell_index_from_ijk(c_x,c_y,c_z);
                Cell *c = dummy_cells[cell_index]->real_cell;

                c->total_pixels++;
                all_pixels++;

                c->pixels += *system->world_grid->get_voxel(i,j,k)<voxel_type_wall;
                filled_pixels += *system->world_grid->get_voxel(i,j,k)<voxel_type_wall;
            }
        }
    }
    porosity = (float)filled_pixels / all_pixels;
}

void ThreadControl::update_mpi(int dimension) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;

    for(int higher_dim=0;higher_dim<=1;higher_dim++) {
        int neighbor_node_index = 2*dimension+higher_dim;
        int node_id = neighbor_nodes[neighbor_node_index];
        if(node_id == myid) continue;

        int num_send = num_molecules_to_be_moved[neighbor_node_index];
        int num_receive = 0;

        if(node_id<myid) {
            MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);
            MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);

            MPI_Send(molecules_to_be_moved[neighbor_node_index],9*num_send,MPI_DOUBLE,node_id,20,MPI_COMM_WORLD);
            MPI_Recv(mpi_receive_buffer,9*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
        }
        else {
            MPI_Recv(&num_receive,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
            MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);

            MPI_Recv(mpi_receive_buffer,9*num_receive,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
            MPI_Send(molecules_to_be_moved[neighbor_node_index],9*num_send,MPI_DOUBLE,node_id,20,MPI_COMM_WORLD);
        }

        memcpy(&new_molecules[9*num_new_molecules], mpi_receive_buffer, 9*num_receive*sizeof(double));
        memset(&molecule_moved[num_new_molecules],0,num_receive*sizeof(bool));
        num_new_molecules += num_receive;

        num_molecules_to_be_moved[neighbor_node_index] = 0;
    }
}

void ThreadControl::add_molecule_to_cell(struct Molecule &molecule, int cell_index) {
    cout << "WARNING, NOT IMPLEMENTED" << endl;
}

void ThreadControl::update_new_molecules(int dimension) {
    for(int n=0;n<num_new_molecules;n++) {
        if(molecule_moved[n]) continue; // Skip molecules that are already moved
        double *r  = &new_molecules[9*n];

        int cell_index = cell_index_from_position(r);
        DummyCell *dummy_cell = dummy_cells[cell_index];

        if(dummy_cell->node_delta_index_vector[dimension] != 0) {
            // We changed cell, and in the dimension we work on right now
            int higher_dim = dummy_cell->node_delta_index_vector[dimension]>0;
            int neighbor_node_id = 2*dimension+higher_dim;
            // Copy the 9 doubles over to the list
            memcpy(&molecules_to_be_moved[neighbor_node_id][ 9*num_molecules_to_be_moved[neighbor_node_id] ],&new_molecules[9*n],9*sizeof(double));
            num_molecules_to_be_moved[neighbor_node_id]++;

            molecule_moved[n] = true;
        }
    }
}

void ThreadControl::add_new_molecules_to_cells() {
    for(int n=0;n<num_new_molecules;n++) {
        if(molecule_moved[n]) continue; // Skip molecules that are already moved
        double *r  = &new_molecules[9*n+0];
        double *v  = &new_molecules[9*n+3];
        double *r0  = &new_molecules[9*n+6];

        int cell_index = cell_index_from_position(r);
        DummyCell *dummy_cell = dummy_cells[cell_index];
        dummy_cell->real_cell->add_molecule(r,v,r0);
    }

    num_new_molecules = 0;
}

ThreadControl::ThreadControl() {

}

