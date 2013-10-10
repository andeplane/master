#include <moviedata.h>
#include <fstream>
#include <iostream>

#ifdef OPENGL
#include <GL/glfw.h>      // Include OpenGL Framework library
#endif

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;

Timestep::Timestep() {
    next = NULL;
    previous = NULL;
}

void Timestep::add_molecule_data(vector<float> &new_positions) {
    positions.insert(positions.end(), new_positions.begin(), new_positions.end());
}

#ifdef OPENGL
void Timestep::render() {
    int num_molecules = positions.size() / 3;
    glColor4f(1.0, 1.0, 1.0, 1.0);
    glPointSize(3.0);
    glBegin(GL_POINTS);
    for(int i=0; i<num_molecules; i++) {
        glVertex3f(10*positions[3*i+0], 10*positions[3*i+1], 10*positions[3*i+2]);
    }
    glEnd();
}
#endif
MovieData::MovieData(int cpus_, int timesteps_)
{
    cpus = cpus_;
    timesteps = timesteps_;
}

void MovieData::load_movie_files(string state_folder) {
    char filename[1000];
    int actual_movie_molecules = 0;
    vector<float> positions;
    Timestep *previous_timestep = NULL;
    for(int timestep=0; timestep < timesteps; timestep++) {
        Timestep *timestep_object = new Timestep();
        if(previous_timestep == NULL) first_timestep = timestep_object;
        else {
            timestep_object->previous = previous_timestep;
            previous_timestep->next = timestep_object;
        }

        previous_timestep = timestep_object;
    }

    // Add periodic time
    first_timestep->previous = previous_timestep;
    previous_timestep->next = first_timestep;

    for(int cpu=0; cpu<cpus; cpu++) {
        sprintf(filename,"%s/movie_files/movie%04d.bin",state_folder.c_str(),cpu);
        ifstream file (filename, ios::in | ios::binary);
        Timestep *timestep_object = first_timestep;
        for(int timestep=0; timestep < timesteps; timestep++) {
            file.read (reinterpret_cast<char*>(&actual_movie_molecules), sizeof(int));
            if(positions.size() < 3*actual_movie_molecules) positions.resize(3*actual_movie_molecules);
            file.read (reinterpret_cast<char*>(&positions[0]), 3*actual_movie_molecules*sizeof(float));
            timestep_object->add_molecule_data(positions);
            timestep_object = timestep_object->next;
        }
        file.close();
    }

    positions.clear();
}
