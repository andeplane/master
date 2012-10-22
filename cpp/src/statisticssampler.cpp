#include "StatisticsSampler.h"
#include <armadillo>
#include "Molecule.h"
#include "Cell.h"
#include <math.h>
#include "Sorter.h"

using namespace arma;

StatisticsSampler::StatisticsSampler(System *_system, int _timesteps) {
    system = _system;
    printTemperature = true;
    temperatureSamples = 0;
    temperatureSum = 0;
    temperatureFile = fopen("temperature.dat","w");

    printVelocityProfile = true;
    velocityProfileSamples = 0;
	
    velocityFile = fopen("velocity.dat","w");
}

void StatisticsSampler::finish() {
	
}

void StatisticsSampler::sample() {
    calculateTemperature();
    calculateVelocityProfile();
}

void StatisticsSampler::calculateTemperature() {
    if(!printTemperature) return;
    temperatureSamples++;

    long N = system->N;
	double energy = 0;
	
    Molecule **molecules = system->molecules;
	for(int n=0;n<N;n++) {

        energy += dot(molecules[n]->v,molecules[n]->v);
	}

    double T = energy/(3*N);

    temperatureSum += T;
	
    fprintf(temperatureFile, "%f %f \n",system->t, T);
}

void StatisticsSampler::calculateVelocityProfile() {
    velocityProfileSamples++;
    if(velocityProfileSamples % 5) return;

	int N = 1000;
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
        fprintf(velocityFile,"%f ",velocities[n]);

    fprintf(velocityFile,"\n");
}











