#include <StatisticsSampler.h>
#include <armadillo>
#include <Molecule.h>
#include <Cell.h>
#include <math.h>
#include <Sorter.h>
#include <CInIFile.h>
using namespace arma;

StatisticsSampler::StatisticsSampler(System *_system, CIniFile &ini) {
    system = _system;
    print_temperature = ini.getbool("print_temperature");
    print_pressure = ini.getbool("print_pressure");
    print_velocity_profile = ini.getbool("print_velocity_profile");
    sample_every_n = ini.getint("sample_every_n");

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
}

void StatisticsSampler::finish() {
	
}

void StatisticsSampler::sample() {
    calculateTemperature();
    calculateVelocityProfile();
}

void StatisticsSampler::calculateTemperature() {
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

void StatisticsSampler::calculatePressure() {
    if(!print_pressure) return;
    // if(++pressure_field_samples % sample_every_n) return;
    pressure_sum = 0;
    representative_cells = 0;
    for(int j=0;j<system->cells_y;j++) {
        for(int i=0;i<system->cells_x;i++) {
            fprintf(pressure_file,"%f ",system->cells[i][j]->pressure);
            if(system->cells[i][j]->pressure > 0) {
                pressure_sum += system->cells[i][j]->pressure;
                representative_cells++;
            }

        }
        fprintf(pressure_file,"\n");
    }
}

void StatisticsSampler::calculateVelocityProfile() {
    if(!print_velocity_profile) return;
    if(++velocity_profile_samples % sample_every_n) return;

    int N = 100;
	Molecule *molecule;
    vec velocities = zeros<vec>(N,1);

	int n;

    for(int i=0;i<system->N;i++) {
        molecule = system->molecules[i];
        n = N*molecule->r(1)/system->height;
        velocities(n) += molecule->v(0)/N;
	}
	
	for(int n=0;n<N;n++) 
        fprintf(velocity_file,"%f ",velocities(n));

    fprintf(velocity_file,"\n");
}
