#include <dsmc_io.h>

#include <system.h>
#include <molecule.h>
#include <cell.h>
#include <mpi.h>
#include <threadcontrol.h>

DSMC_IO::DSMC_IO(System *system_) {
    system = system_;
    settings = system->settings;
    movie_frames = 0;
    movie_file_open = false;
    if(system->myid==0) energy_file = fopen("energy.txt","w");
}

void DSMC_IO::save_state_to_movie_file() {
    if(settings->create_movie && !(system->steps % settings->movie_every_n_frame)) {
        if(!movie_file_open) {
            char *filename = new char[100];
            sprintf(filename,"state_files/movie%04d.bin",system->myid);
            movie_file = new ofstream(filename,ios::out | ios::binary);
            movie_file_open = true;
            data = new double[3*system->thread_control.allocated_particle_data];
            delete filename;
        }

        int count = 0;
        for(unsigned int i=0;i<system->thread_control.cells.size();i++) {
            Cell *c = system->thread_control.cells[i];
            for(unsigned int j=0;j<c->molecules.size();j++) {
                Molecule *m = c->molecules[j];
                data[count++] = m->r[0];
                data[count++] = m->r[1];
                data[count++] = m->r[2];
            }
        }

        count /= 3; // This should represent the number of particles

        movie_file->write (reinterpret_cast<char*>(&count), sizeof(int));
        movie_file->write (reinterpret_cast<char*>(data), 3*count*sizeof(double));
    }
}

void DSMC_IO::save_state_to_file_binary() {
    if(system->myid==0) cout << "Saving state to file..." << endl;
    ThreadControl &thread_control = system->thread_control;

    int N = thread_control.num_particles;

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);

    ofstream file (filename, ios::out | ios::binary);
    double *tmp_data = new double[9*N];

    int count = 0;
    for(unsigned int i=0;i<thread_control.cells.size();i++) {
        Cell *c = thread_control.cells[i];
        for(unsigned int j=0;j<c->molecules.size();j++) {
            Molecule *m = c->molecules[j];
            tmp_data[count++] = m->r[0];
            tmp_data[count++] = m->r[1];
            tmp_data[count++] = m->r[2];

            tmp_data[count++] = m->r_initial[0];
            tmp_data[count++] = m->r_initial[1];
            tmp_data[count++] = m->r_initial[2];

            tmp_data[count++] = m->v[0];
            tmp_data[count++] = m->v[1];
            tmp_data[count++] = m->v[2];
        }
    }

    file.write (reinterpret_cast<char*>(&N), sizeof(int));
    file.write (reinterpret_cast<char*>(tmp_data), 9*N*sizeof(double));

    file.close();
    delete tmp_data;
    delete filename;
}

void DSMC_IO::load_state_from_file_binary() {
    if(system->myid==0) cout << "Loading state from file..." << endl;
    int N = 0;
    ThreadControl &thread_control = system->thread_control;

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);
    ifstream file (filename, ios::in | ios::binary);

    file.read(reinterpret_cast<char*>(&N),sizeof(int));

    double *tmp_data = new double[9*N];

    file.read(reinterpret_cast<char*>(tmp_data), 9*N*sizeof(double));
    file.close();

    for(int n=0;n<N;n++) {
        Molecule *m = new Molecule(system);
        m->atoms = system->eff_num;
        m->r = &thread_control.positions[3*thread_control.all_molecules.size()];
        m->v = &thread_control.velocities[3*thread_control.all_molecules.size()];
        m->r_initial = &thread_control.initial_positions[3*thread_control.all_molecules.size()];
        m->r[0] = tmp_data[9*n+0];
        m->r[1] = tmp_data[9*n+1];
        m->r[2] = tmp_data[9*n+2];

        m->r_initial[0] = tmp_data[9*n+3];
        m->r_initial[1] = tmp_data[9*n+4];
        m->r_initial[2] = tmp_data[9*n+5];

        m->v[0] = tmp_data[9*n+6];
        m->v[1] = tmp_data[9*n+7];
        m->v[2] = tmp_data[9*n+8];

        thread_control.dummy_cells[thread_control.cell_index_from_molecule(m)]->real_cell->add_molecule(m);
        thread_control.all_molecules.push_back(m);
    }

    delete filename;
    delete tmp_data;
}

void DSMC_IO::finalize() {
    if(movie_file_open) {
        movie_file->close();
    }

    if(system->myid != 0) return;
    fclose(energy_file);
}
