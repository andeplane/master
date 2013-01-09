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
    print_pressure = ini->getbool("print_pressure");
    print_velocity_profile = ini->getbool("print_velocity_profile");
    sample_every_n = ini->getint("sample_every_n");
    print_velocity_field = ini->getbool("print_velocity_field");

    temperature_samples = 0;
    temperature_sum = 0;
    if(print_temperature)
        temperature_file = fopen("temperature.dat","w");

    pressure_samples = 0;
    pressure_sum = 0;
    if(print_pressure)
        pressure_file = fopen("pressure.dat","w");

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

void StatisticsSampler::calculate_pressure() {
    if(!print_pressure) return;
    // if(++pressure_field_samples % sample_every_n) return;
    pressure_sum = 0;
    representative_cells = 0;
    for(int j=0;j<system->cells_y;j++) {
        for(int i=0;i<system->cells_x;i++) {
            fprintf(pressure_file,"%f ",system->cells[i][j]->pressure);
            if(system->cells[i][j]->pressure > 0) {
                pressure_sum += system->cells[i][j]->f_sum;
                representative_cells++;
            }
        }
        fprintf(pressure_file,"\n");
    }
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
double StatisticsSampler::get_pressure() {
    if(!print_pressure) return 0;
    double P = (system->N*system->T + 1.0/2*pressure_sum/representative_cells)/system->volume;

    return P;
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
