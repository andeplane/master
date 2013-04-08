#include <dsmc_io.h>

#include <system.h>
#include <cell.h>
#include <mpi.h>
#include <threadcontrol.h>
#include <grid.h>
#include <settings.h>
DSMC_IO::DSMC_IO(System *system_) {
    system = system_;
    settings = system->settings;
    movie_frames = 0;
    movie_file_open = false;
    if(system->myid==0) {
        energy_file = fopen("energy.txt","w");
        velocity_file = fopen("velocity.txt","w");
    }
}

void DSMC_IO::save_state_to_movie_file() {
    if(settings->create_movie && !(system->steps % settings->movie_every_n_frame)) {
        if(!movie_file_open) {
            char *filename = new char[100];
            sprintf(filename,"state_files/movie%04d.bin",system->myid);
            movie_file = new ofstream(filename,ios::out | ios::binary);
            movie_file_open = true;
            data = new double[3*MAX_MOLECULE_NUM];
            delete filename;
        }

        int count = 0;
        for(unsigned int n=0;n<system->thread_control.num_molecules;n++) {
            data[count++] = system->thread_control.r[3*n+0];
            data[count++] = system->thread_control.r[3*n+1];
            data[count++] = system->thread_control.r[3*n+2];
        }

        count /= 3; // This should represent the number of particles

        movie_file->write (reinterpret_cast<char*>(&count), sizeof(int));
        movie_file->write (reinterpret_cast<char*>(data), 3*count*sizeof(double));
    }
}

void DSMC_IO::save_state_to_file_binary() {
    if(system->myid==0) cout << "Saving state to file..." << endl;
    ThreadControl &thread_control = system->thread_control;

    int N = thread_control.num_molecules;

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);

    ofstream file (filename, ios::out | ios::binary);

    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << endl;
        exit(1);
    }

    double *tmp_data = new double[9*N];

    int count = 0;
    for(unsigned int n=0;n<thread_control.num_molecules;n++) {

        tmp_data[count++] = thread_control.r[3*n+0];
        tmp_data[count++] = thread_control.r[3*n+1];
        tmp_data[count++] = thread_control.r[3*n+2];

        tmp_data[count++] = thread_control.v[3*n+0];
        tmp_data[count++] = thread_control.v[3*n+1];
        tmp_data[count++] = thread_control.v[3*n+2];

        tmp_data[count++] = thread_control.r0[3*n+0];
        tmp_data[count++] = thread_control.r0[3*n+1];
        tmp_data[count++] = thread_control.r0[3*n+2];
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
    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << endl;
        exit(1);
    }

    file.read(reinterpret_cast<char*>(&N),sizeof(int));
    double *tmp_data = new double[9*N];

    file.read(reinterpret_cast<char*>(tmp_data), 9*N*sizeof(double));
    file.close();

    for(int n=0;n<N;n++) {
        double *r = &thread_control.r[3*n];
        double *v = &thread_control.v[3*n];
        double *r0 = &thread_control.r0[3*n];

        r[0] = tmp_data[9*n+0];
        r[1] = tmp_data[9*n+1];
        r[2] = tmp_data[9*n+2];
        v[0] = tmp_data[9*n+3];
        v[1] = tmp_data[9*n+4];
        v[2] = tmp_data[9*n+5];
        r0[0] = tmp_data[9*n+6];
        r0[1] = tmp_data[9*n+7];
        r0[2] = tmp_data[9*n+8];

        Cell *cell = thread_control.all_cells[thread_control.cell_index_from_position(r)];
        cell->add_molecule(n,thread_control.molecule_index_in_cell,thread_control.molecule_cell_index);
    }

    thread_control.num_molecules = N;

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

void DSMC_IO::read_grid_matrix(string filename, Grid *grid) {
    ifstream file (filename.c_str(), ios::in | ios::binary);
    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << endl;
        exit(1);
    }
    int Nx, Ny, Nz, points;

    file.read (reinterpret_cast<char*>(&Nx), sizeof(int));
    file.read (reinterpret_cast<char*>(&Ny), sizeof(int));
    file.read (reinterpret_cast<char*>(&Nz), sizeof(int));
    points = Nx*Ny*Nz;
    grid->Nx = Nx; grid->Ny = Ny; grid->Nz = Nz; grid->points = points;

    grid->voxels = new unsigned char[points];
    grid->normals   = new float[3*points];
    grid->tangents1 = new float[3*points];
    grid->tangents2 = new float[3*points];

    file.read (reinterpret_cast<char*>(grid->voxels), points*sizeof(unsigned char));
    file.read (reinterpret_cast<char*>(grid->normals), 3*points*sizeof(float));
    file.read (reinterpret_cast<char*>(grid->tangents1), 3*points*sizeof(float));
    file.read (reinterpret_cast<char*>(grid->tangents2), 3*points*sizeof(float));
    file.close();
}
