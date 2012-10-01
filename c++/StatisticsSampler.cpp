#include "StatisticsSampler.h"
#include <armadillo>
#include "Molecule.h"
#include "Cell.h"
#include <math.h>
#include "Sorter.h"

using namespace arma;

StatisticsSampler::StatisticsSampler(System *system, int timesteps) {
	this->system = system;
	this->printTemperature = true;
	this->temperatureSamples = 0;
	this->temperatureSum = 0;
	this->temperatureFile = fopen("temperature.dat","w");
}

void StatisticsSampler::sample() {
	this->calculateTemperature();
}

void StatisticsSampler::calculateTemperature() {
	if(!this->printTemperature) return;
	this->temperatureSamples++;

	long N = this->system->N;
	double energy = 0;
	
	Molecule **molecules = this->system->molecules;
	for(int n=0;n<N;n++) {
		energy += dot(molecules[n]->v,molecules[n]->v);
	}

	double T = energy/(3*N);

	this->temperatureSum += T;
	
	fprintf(this->temperatureFile, "%f %f \n",this->system->t, T);
}
