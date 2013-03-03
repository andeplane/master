#include <dsmc_io.h>
#include <fstream>

using namespace std;

DSMC_IO::DSMC_IO(Settings *settings_, System *system_) {
    settings = settings_;
    system = system_;
}


void DSMC_IO::state_to_file_binary() {
    time_t t0;
    t0 = clock();
    cout << "Saving state to file..." << endl;

    int N = system->N;
    ofstream file ("state", ios::out | ios::binary);

    file.write (reinterpret_cast<char*>(&N), sizeof(int));
    file.write (reinterpret_cast<char*>(system->positions), 3*N*sizeof(double));
    file.write (reinterpret_cast<char*>(system->initial_positions), 3*N*sizeof(double));
    file.write (reinterpret_cast<char*>(system->velocities), 3*N*sizeof(double));

    file.close();

    double t = ((double)clock()-t0)/CLOCKS_PER_SEC;
    cout << "Saving took " << t << " seconds." << endl;
}

void DSMC_IO::state_from_file_binary() {
    time_t t0;
    t0 = clock();
    cout << "Loading state from file..." << endl;
    int N = 0;
    ifstream file ("state", ios::in | ios::binary);
    file.read(reinterpret_cast<char*>(&N),sizeof(int));

    system->N = N;
    system->positions = new double[3*N];
    system->velocities = new double[3*N];
    system->initial_positions = new double[3*N];
    system->molecules.reserve(N);

    file.read(reinterpret_cast<char*>(system->positions), 3*N*sizeof(double));
    file.read(reinterpret_cast<char*>(system->initial_positions), 3*N*sizeof(double));
    file.read(reinterpret_cast<char*>(system->velocities), 3*N*sizeof(double));
    file.close();

    for(int n=0;n<N;n++) {
        Molecule *m = new Molecule(system);
        m->atoms = system->eff_num;
        m->index = n;
        m->r = &system->positions[3*n];
        m->v = &system->velocities[3*n];
        m->initial_r = &system->initial_positions[3*n];
        system->molecules.push_back(m);
    }

    double t = ((double)clock()-t0)/CLOCKS_PER_SEC;
    cout << "Loading took " << t << " seconds." << endl;
}
