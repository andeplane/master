#include <mdio.h>
#include <system.h>
#include <settings.h>
#include <mdtimer.h>

MDIO::MDIO()
{

}

void MDIO::setup(System *system_) {
    system = system_;
    settings = system->settings;
    movie_file_open = false;
    if(system->myid==0) {
        energy_file = fopen("statistics/energy.txt","w");
        pressure_file = fopen("statistics/pressure.txt","w");
        velocity_file = fopen("statistics/velocity.txt","w");
    }
}

void MDIO::save_state_to_movie_file() {
    system->mdtimer->start_io();
    if(settings->create_movie && !(system->steps % settings->movie_every_n_frame)) {
        if(!movie_file_open) {
            char *filename = new char[100];
            sprintf(filename,"movie_files/movie%04d.bin",system->myid);
            movie_file = new ofstream(filename,ios::out | ios::binary);
            movie_file_open = true;
            data = new double[3*MAX_ATOM_NUM];

            delete filename;
        }

        for(unsigned long n=0;n<system->num_atoms_local;n++) {
            data[3*n+0] = system->positions[n][0] + system->origo[0];
            data[3*n+1] = system->positions[n][1] + system->origo[1];
            data[3*n+2] = system->positions[n][2] + system->origo[2];
        }

        movie_file->write (reinterpret_cast<char*>(&system->num_atoms_local), sizeof(unsigned long));
        movie_file->write (reinterpret_cast<char*>(data), 3*system->num_atoms_local*sizeof(double));
        movie_file->write (reinterpret_cast<char*>(&system->atom_type), system->num_atoms_local*sizeof(unsigned long));
        movie_file->flush();
    }
    system->mdtimer->end_io();
}

void MDIO::save_state_to_file_binary() {
    system->mdtimer->start_io();

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);

    ofstream file (filename, ios::out | ios::binary);
    double *tmp_data = new double[6*system->num_atoms_local];

    for(unsigned int i=0;i<system->num_atoms_local;i++) {
        tmp_data[6*i + 0] = system->positions[i][0] + system->origo[0];
        tmp_data[6*i + 1] = system->positions[i][1] + system->origo[1];
        tmp_data[6*i + 2] = system->positions[i][2] + system->origo[2];

        tmp_data[6*i + 3] = system->velocities[3*i+0];
        tmp_data[6*i + 4] = system->velocities[3*i+1];
        tmp_data[6*i + 5] = system->velocities[3*i+2];
    }

    file.write (reinterpret_cast<char*>(&system->num_atoms_local), sizeof(int));
    file.write (reinterpret_cast<char*>(tmp_data), 6*system->num_atoms_local*sizeof(double));
    file.write (reinterpret_cast<char*>(&system->atom_type), system->num_atoms_local*sizeof(unsigned long));

    file.close();
    delete tmp_data;
    delete filename;
    system->mdtimer->end_io();
}

void MDIO::load_state_from_file_binary() {
    system->mdtimer->start_io();

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);

    ifstream file (filename, ios::in | ios::binary);

    file.read(reinterpret_cast<char*>(&system->num_atoms_local),sizeof(int));

    double *tmp_data = new double[6*system->num_atoms_local];
    file.read(reinterpret_cast<char*>(tmp_data),6*system->num_atoms_local*sizeof(double));
    file.read(reinterpret_cast<char*>(&system->atom_type), system->num_atoms_local*sizeof(unsigned long));
    file.close();

    for(unsigned int i=0;i<system->num_atoms_local;i++) {
        system->positions[i][0] = tmp_data[6*i+0] - system->origo[0];
        system->positions[i][1] = tmp_data[6*i+1] - system->origo[1];
        system->positions[i][2] = tmp_data[6*i+2] - system->origo[2];

        system->velocities[3*i+0] = tmp_data[6*i+3];
        system->velocities[3*i+1] = tmp_data[6*i+4];
        system->velocities[3*i+2] = tmp_data[6*i+5];
    }
    delete tmp_data;
    delete filename;
    system->mdtimer->end_io();
}

void MDIO::finalize() {
    if(movie_file_open) {
        movie_file->close();
    }

    if(system->myid != 0) return;
    fclose(energy_file);
    fclose(pressure_file);
    fclose(velocity_file);
}
