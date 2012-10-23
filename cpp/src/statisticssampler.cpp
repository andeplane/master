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
    print_velocity_profile = ini.getbool("print_velocity_profile");

    temperature_samples = 0;
    temperature_sum = 0;

    temperature_file = fopen("temperature.dat","w");

    velocity_profile_samples = 0;
	
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
    if(++temperature_samples % 5) return;

    long N = system->N;
	double energy = 0;
	
    Molecule **molecules = system->molecules;
	for(int n=0;n<N;n++) {
        energy += dot(molecules[n]->v,molecules[n]->v);
	}

    double T = energy/(3*N);

    temperature_sum += T;
	
    fprintf(temperature_file, "%f %f \n",system->t, T);
}

void StatisticsSampler::calculateVelocityProfile() {
    if(!print_velocity_profile) return;
    if(++velocity_profile_samples % 5) return;

    int N = 100;
	Molecule *molecule;
	double *velocities = new double[100];
	for(int n=0;n<N;n++) 
		velocities[n] = 0;

	int n;
	for(int i=0;i<this->system->N;i++) {
        molecule = system->molecules[i];
        n = N*molecule->r(1)/system->height;
        velocities[n] += molecule->v(0);
	}
	
	for(int n=0;n<N;n++) 
        fprintf(velocity_file,"%f ",velocities[n]);

    fprintf(velocity_file,"\n");
}
