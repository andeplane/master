#include <threadcontrol.h>
#include <system.h>
#include <settings.h>
#include <cell.h>
#include <molecule.h>
#include <random.h>
#include <grid.h>
#include <mpi.h>

void ThreadControl::setup(System *system_) {
    system = system_;
    settings = system->settings;
    num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;
    myid = system->myid;

    idx = myid/(settings->nodes_y*settings->nodes_z); // Node id in x-direction
    idy = (myid/settings->nodes_z) % settings->nodes_y; // Node id in y-direction
    idz = myid % settings->nodes_z; // Node id in z-direction

    origo(0) = idx*(system->Lx/settings->nodes_x);
    origo(1) = idy*(system->Ly/settings->nodes_y);
    origo(2) = idz*(system->Lz/settings->nodes_z);
    if(myid==0) cout << "Setting up cells..." << endl;
    setup_cells();
    if(myid==0) cout << "Updating cell volume..." << endl;
    update_cell_volume();
    if(myid==0) cout << "Setting up molecules..." << endl;
    setup_molecules();
    mpi_data = new double[9*2*settings->number_of_particles];
}

void ThreadControl::update_cell_volume() {
    for(int i=0;i<cells.size();i++) {
        Cell *c = cells[i];
        c->vr_max = 3*system->mpv;
        c->update_volume();
    }
}

void ThreadControl::setup_molecules() {
    free_molecules.reserve(1000);
    num_particles = settings->number_of_particles*porosity/num_nodes;

    for(int n=0;n<num_particles;n++) {
        Molecule *m = new Molecule(system);
        m->atoms = system->eff_num;
        m->v(0) = 0.6; //system->rnd->nextGauss()*sqrt(system->temperature);
        m->v(1) = 0; //system->rnd->nextGauss()*sqrt(system->temperature);
        m->v(2) = 0; //system->rnd->nextGauss()*sqrt(system->temperature);
        find_position(m);
        dummy_cells[cell_index_from_molecule(m)]->real_cell->add_molecule(m);
    }
}

inline void ThreadControl::find_position(Molecule *m) {
    bool did_collide = true;
    while(did_collide) {
        m->r(0) = origo(0) + system->Lx/settings->nodes_x*system->rnd->nextDouble();
        m->r(1) = origo(1) + system->Ly/settings->nodes_y*system->rnd->nextDouble();
        m->r(2) = origo(2) + system->Lz/settings->nodes_z*system->rnd->nextDouble();

        did_collide = *system->world_grid->get_voxel(m->r)>=voxel_type_wall;
    }

    m->r_initial(0) = m->r(0);
    m->r_initial(1) = m->r(1);
    m->r_initial(2) = m->r(2);
}

inline int ThreadControl::cell_index_from_ijk(const int &i, const int &j, const int &k) {
    return i*system->cells_y*system->cells_z + j*system->cells_z + k;
}

int ThreadControl::cell_index_from_molecule(Molecule *m) {
    int i = m->r(0)/system->Lx*system->cells_x;
    int j = m->r(1)/system->Ly*system->cells_y;
    int k = m->r(2)/system->Lz*system->cells_z;

    return cell_index_from_ijk(i,j,k);
}

void ThreadControl::setup_cells() {
    int num_cells_global = num_nodes*settings->cells_per_node_x*settings->cells_per_node_y*settings->cells_per_node_z;
    nodes_new_atoms_list.resize(num_nodes);
    dummy_cells.reserve(num_cells_global);

    for(int i=0;i<system->cells_x;i++) {
        for(int j=0;j<system->cells_y;j++) {
            for(int k=0;k<system->cells_z;k++) {
                int c_idx = i/settings->cells_per_node_x;
                int c_idy = j/settings->cells_per_node_y;
                int c_idz = k/settings->cells_per_node_z;
                int node_id = c_idx*settings->nodes_y*settings->nodes_z + c_idy*settings->nodes_z + c_idz;

                DummyCell *dc = new DummyCell();
                dc->node_id = node_id;
                dc->index = cell_index_from_ijk(i,j,k);
                // dc->new_molecules.reserve(1000);
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

    int i_start = float(this->idx)/system->cells_x*system->world_grid->Nx;
    int i_end   = float(this->idx+1)/system->cells_x*system->world_grid->Nx;

    int j_start = float(this->idy)/system->cells_y*system->world_grid->Ny;
    int j_end   = float(this->idy+1)/system->cells_y*system->world_grid->Ny;

    int k_start = float(this->idz)/system->cells_z*system->world_grid->Nz;
    int k_end   = float(this->idz+1)/system->cells_z*system->world_grid->Nz;
    int cell_index, c_x, c_y, c_z;

    for(int k=k_start;k<k_end;k++) {
        for(int j=j_start;j<j_end;j++) {
            for(int i=i_start;i<i_end;i++) {
                c_x = (i-i_start)*settings->cells_per_node_x/(float)system->world_grid->Nx;
                c_y = (j-j_start)*settings->cells_per_node_y/(float)system->world_grid->Ny;
                c_z = (k-k_start)*settings->cells_per_node_z/(float)system->world_grid->Nz;

                cell_index = c_x + c_y*settings->cells_per_node_x + c_z*settings->cells_per_node_x*settings->cells_per_node_y;
                Cell *c = cells[cell_index];

                c->total_pixels++;
                all_pixels++;

                c->pixels += *system->world_grid->get_voxel(i,j,k)<voxel_type_wall;
                filled_pixels += *system->world_grid->get_voxel(i,j,k)<voxel_type_wall;
            }
        }
    }
    porosity = (float)filled_pixels / all_pixels;
}

void ThreadControl::update_local_cells() {
    for(int i=0;i<cells.size();i++) {
        Cell *c = cells[i];
        DummyCell *dc = c->dummy_cell;
        for(int j=0;j<dc->new_molecules.size();j++) {
            Molecule *m = dc->new_molecules[j];
            // Remove from old cell
            dummy_cells[m->cell_index]->real_cell->remove_molecule(m);
            c->add_molecule(m);
        }
        dc->new_molecules.clear();
    }
}

double comm_time = 0;

void ThreadControl::update_mpi() {
    MPI_Barrier(MPI_COMM_WORLD);
    double t0;
    if(myid==0) t0=MPI_Wtime();

    // cout << "Mpi update staring on node " << myid << endl;

    for(int i=0;i<num_nodes;i++) {
        vector<Molecule*> &molecules = nodes_new_atoms_list[i];

        if(i==myid) {
            for(int j=0;j<num_nodes;j++) {
                if(j==myid) continue;

                int num_received;
                MPI_Recv(&num_received, 1, MPI_INT, j, 100,
                             MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                if(num_received) {
                    MPI_Recv(mpi_data, num_received, MPI_DOUBLE, j, 100,
                                 MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    for(int n=0;n<num_received/9;n++) {
                        Molecule *m;
                        if(free_molecules.size()>0) {
                            m = free_molecules.back();
                            free_molecules.erase(free_molecules.end()-1);
                        } else {
                            m = new Molecule(system);
                        }

                        m->atoms = system->eff_num;
                        m->r(0) = mpi_data[9*n+0];
                        m->r(1) = mpi_data[9*n+1];
                        m->r(2) = mpi_data[9*n+2];

                        m->v(0) = mpi_data[9*n+3];
                        m->v(1) = mpi_data[9*n+4];
                        m->v(2) = mpi_data[9*n+5];

                        m->r_initial(0) = mpi_data[9*n+6];
                        m->r_initial(1) = mpi_data[9*n+7];
                        m->r_initial(2) = mpi_data[9*n+8];

                        dummy_cells[cell_index_from_molecule(m)]->real_cell->add_molecule(m);
                        num_particles++;
                    }
                }
            }

            continue;
        }
        int count = 0;

        for(int n=0;n<molecules.size();n++) {
            Molecule *m = molecules[n];
            mpi_data[count++] = m->r(0);
            mpi_data[count++] = m->r(1);
            mpi_data[count++] = m->r(2);
            mpi_data[count++] = m->v(0);
            mpi_data[count++] = m->v(1);
            mpi_data[count++] = m->v(2);
            mpi_data[count++] = m->r_initial(0);
            mpi_data[count++] = m->r_initial(1);
            mpi_data[count++] = m->r_initial(2);

            dummy_cells[m->cell_index]->real_cell->remove_molecule(m);
            num_particles--;
            free_molecules.push_back(m);
        }

        molecules.clear();

        MPI_Send(&count, 1, MPI_INT, i, 100,
                     MPI_COMM_WORLD);
        if(count) {
            MPI_Send(mpi_data, count, MPI_DOUBLE, i, 100,
                         MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid==0) comm_time += MPI_Wtime()-t0;
    if(!(system->steps%500)) {
        if(myid==0) cout << "MPI-COMM TOTAL: " << comm_time << " seconds." << endl;
    }
}

ThreadControl::ThreadControl() {

}
