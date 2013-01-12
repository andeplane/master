#include <statisticssampler.h>
#include <armadillo>
#include <molecule.h>
#include <cell.h>
#include <math.h>
#include <sorter.h>
#include <CIniFile.h>
using namespace arma;

StatisticsSampler::StatisticsSampler(System *_system, CIniFile *_ini) {
    system = _system;
    ini = _ini;
    print_temperature = ini->getbool("print_temperature");
    print_velocity_profile = ini->getbool("print_velocity_profile");
    sample_every_n = ini->getint("sample_every_n");
    print_velocity_field = ini->getbool("print_velocity_field");

    temperature_samples = 0;
    temperature_sum = 0;
    if(print_temperature)
        temperature_file = fopen("temperature.dat","w");

    velocity_profile_samples = 0;
    if(print_velocity_profile)
        velocity_file = fopen("velocity.dat","w");

    if(print_velocity_field) {
        velocity_field_file_x = fopen("vel_x.dat","w");
        velocity_field_file_y = fopen("vel_y.dat","w");;
    }
}

void StatisticsSampler::finish() {
	
}

void StatisticsSampler::sample() {
    calculate_temperature();
    calculate_velocity_profile();
    update_cell_statistics();
}

void StatisticsSampler::calculate_temperature() {
    if(!print_temperature) return;
    if(++temperature_samples % sample_every_n) return;

    long N = system->N;
    double T = 0;
	
    Molecule **molecules = system->molecules;
	for(int n=0;n<N;n++) {
        T += dot(molecules[n]->v,molecules[n]->v)/(3*N);
	}
    // Note that E=0.5*M*v^2 and M=N*m, but it cancels in the denominator

    temperature_sum += T;
	
    fprintf(temperature_file, "%f %f \n",system->t, T);
}

void StatisticsSampler::calculate_velocity_profile() {
    if(!print_velocity_profile) return;
    if(++velocity_profile_samples % sample_every_n) return;

    int N = 100;
	Molecule *molecule;
    vec velocities = zeros<vec>(N,1);
    vec velocity_count= zeros<vec>(N,1);

	int n;


    for(int i=0;i<system->N;i++) {
        molecule = system->molecules[i];
        n = N*molecule->r(1)/system->height;
        velocities(n) += molecule->v(0);
        velocity_count(n)++;
	}

    for(n=0;n<N;n++) {
        if(velocity_count(n)>0)
            velocities(n) /= velocity_count(n);
    }
	
	for(int n=0;n<N;n++) 
        fprintf(velocity_file,"%f ",velocities(n));

    fprintf(velocity_file,"\n");
}

double StatisticsSampler::get_temperature() {
    if(!print_temperature) return 0;

    return temperature_sum/(temperature_samples/ini->getint("sample_every_n"));
}

void StatisticsSampler::calculate_velocity_field() {
    if(!print_velocity_field) return;

    Cell *c;
    for(int j=0;j<system->cells_y;j++) {
        for(int i=0;i<system->cells_x;i++) {
            c = system->cells[i][j];

            fprintf(velocity_field_file_x,"%.5f ",c->momentum(0));
            fprintf(velocity_field_file_y,"%.5f ",c->momentum(1));
        }

        fprintf(velocity_field_file_x,"\n");
        fprintf(velocity_field_file_y,"\n");
    }
}

double StatisticsSampler::calculate_diffusion_constant() {
    Molecule *molecule;
    double r_squared = 0;
    for(int i=0;i<system->N;i++) {
        molecule = system->molecules[i];
        r_squared += molecule->squared_distance_from_initial_position();
    }
    r_squared /= system->N;

    double diffusion_constant = r_squared/6/system->t;
    return diffusion_constant;
}


mat StatisticsSampler::calculate_global_pressure_tensor() {
    vector<Molecule*> molecules;
    for(int n=0;n<system->N;n++)
        molecules.push_back(system->molecules[n]);

    return calculate_pressure_tensor(molecules);
}

vector<mat> StatisticsSampler::calculate_local_pressure_tensor() {
    vector<mat> pressure_tensors;

    Cell *c;
    for(int i=0;i<system->cells_x;i++) {
        for(int j=0;j<system->cells_y;j++) {
            c = system->cells[i][j];
            vector<Molecule*> molecules;
            for(int n=0;n<c->particles;n++)
                molecules.push_back(system->molecules[c->particle_indices[n]]);
            pressure_tensors.push_back(calculate_pressure_tensor(molecules));
        }
    }

    return pressure_tensors;
}

mat StatisticsSampler::calculate_pressure_tensor(vector<Molecule*> molecules) {
    mat p = zeros<mat>(3,3);
    vec v = zeros<vec>(3,1);
    double rho = 0;
    int N = molecules.size();
    if(!N) return p;

    Molecule *m;

    for(int n=0;n<molecules.size();n++) {
        m = molecules[n];
        v += m->mass*m->v/N;
        rho += m->mass/N;
        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                p(i,j) += m->v(i)*m->v(j)*m->mass/N;
            }
        }
    }
    v /= rho;

    for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            p(i,j) -= v(i)*v(j)*rho;
        }
    }

    return p;
}

void StatisticsSampler::update_cell_statistics() {
    Cell *c; Molecule *m;
    for(int i=0;i<system->cells_x;i++) {
        for(int j=0;j<system->cells_y;j++) {
            c = system->cells[i][j];
            vector<Molecule*> molecules;
            vec momentum = zeros<vec>(3,1);
            double energy = 0;

            for(int n=0;n<c->particles;n++) {
                m = system->molecules[c->particle_indices[n]];

                energy += 0.5*dot(m->v,m->v)*m->mass*m->atoms;
                momentum += m->v*m->mass*m->atoms;

                molecules.push_back(m);
            }
            c->update_pressure_tensor(calculate_pressure_tensor(molecules));
            c->update_momentum(momentum);
            c->update_energy(energy);
            c->update_temperature(2*energy/(3*c->particles));
        }
    }
}
