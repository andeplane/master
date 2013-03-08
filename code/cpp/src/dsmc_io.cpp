#include <dsmc_io.h>

#include <fstream>
#include <system.h>
#include <molecule.h>
#include <cell.h>

using namespace std;

DSMC_IO::DSMC_IO(System *system_) {
    system = system_;
    settings = system->settings;
    movie_frames = 0;
    movie_file_open = false;
    energy_file = fopen("energy.txt","w");
}

void DSMC_IO::save_state_to_file_xyz() {
//    time_t t0;
//    t0 = clock();
//    cout << "Saving state to xyz-file..." << endl;

//    ofstream file ("state.xyz", ios::out);
//    file << system->N << endl;
//    file << "sup" << endl;

//    for(int n=0;n<system->N;n++) {
//        // We return height - r(1) because system is inverted
//        file << "H " << system->molecules[n]->r[0] << " " << system->molecules[n]->r[1] << " " << system->molecules[n]->r[2] << endl;
//    }

//    file.close();

//    double t = ((double)clock()-t0)/CLOCKS_PER_SEC;
//    cout << "Saving took " << t << " seconds." << endl;
}

void DSMC_IO::save_state_to_movie_file() {
    if(settings->create_movie && !(system->steps % settings->movie_every_n_frame)) {
        if(!movie_file_open) {
            movie_file = fopen("movie.xyz","w");
            movie_file_open = true;
        }

        fprintf(movie_file,"%d\n",system->thread_control.num_particles);
        fprintf(movie_file,"Random comment that must be here\n");

        for(int i=0;i<system->thread_control.cells.size();i++) {
            Cell *c = system->thread_control.cells[i];
            for(int j=0;j<c->molecules.size();j++) {
                Molecule *m = c->molecules[j];
                fprintf(movie_file,"H %.10f %.10f %.10f\n", m->r(0),m->r(1), m->r(2));
            }
        }
    }
}

void DSMC_IO::save_state_to_file_binary() {
//    time_t t0;
//    t0 = clock();
//    cout << "Saving state to file..." << endl;

//    int N = system->N;
//    ofstream file ("state", ios::out | ios::binary);

//    file.write (reinterpret_cast<char*>(&N), sizeof(int));
//    file.write (reinterpret_cast<char*>(system->positions), 3*N*sizeof(double));
//    file.write (reinterpret_cast<char*>(system->initial_positions), 3*N*sizeof(double));
//    file.write (reinterpret_cast<char*>(system->velocities), 3*N*sizeof(double));

//    file.close();

//    double t = ((double)clock()-t0)/CLOCKS_PER_SEC;
//    cout << "Saving took " << t << " seconds." << endl;
}

void DSMC_IO::load_state_from_file_binary() {
//    time_t t0;
//    t0 = clock();
//    cout << "Loading state from file..." << endl;
//    int N = 0;
//    ifstream file ("state", ios::in | ios::binary);
//    file.read(reinterpret_cast<char*>(&N),sizeof(int));

//    system->N = N;
//    system->positions = new double[3*N];
//    system->velocities = new double[3*N];
//    system->initial_positions = new double[3*N];
//    system->molecules.reserve(N);

//    file.read(reinterpret_cast<char*>(system->positions), 3*N*sizeof(double));
//    file.read(reinterpret_cast<char*>(system->initial_positions), 3*N*sizeof(double));
//    file.read(reinterpret_cast<char*>(system->velocities), 3*N*sizeof(double));
//    file.close();

//    for(int n=0;n<N;n++) {
//        Molecule *m = new Molecule(system);
//        m->atoms = system->eff_num;
//        m->r = &system->positions[3*n];
//        m->v = &system->velocities[3*n];
//        m->initial_r = &system->initial_positions[3*n];
//        system->molecules.push_back(m);
//    }

//    double t = ((double)clock()-t0)/CLOCKS_PER_SEC;
//    cout << "Loading took " << t << " seconds." << endl;
}

void DSMC_IO::finalize() {
    if(movie_file_open) {
        fclose(movie_file);
    }
    fclose(energy_file);
}
