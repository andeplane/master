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

	this->printVelocityProfile = true;
	this->velocityProfileSamples = 0;
	
	this->velocityFile = fopen("velocity.dat","w");
}

void StatisticsSampler::finish() {
	
}

void StatisticsSampler::sample() {
	this->calculateTemperature();
	this->calculateVelocityProfile();
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

inline int calcCellIndex(int i, int j, int k, int cellsPerDimension) {
	return k*cellsPerDimension*cellsPerDimension + j*cellsPerDimension + i;
}

void StatisticsSampler::calculateVelocityProfile() {
	this->velocityProfileSamples++;
	if(this->velocityProfileSamples % 5) return;

	int N = 1000;
	Molecule *molecule;
	double *velocities = new double[100];
	for(int n=0;n<N;n++) 
		velocities[n] = 0;

	int n;
	for(int i=0;i<this->system->N;i++) {
		molecule = this->system->molecules[i];
		n = N*molecule->r(1)/this->system->L;
		velocities[n] += molecule->v(0);
	}
	
	for(int n=0;n<N;n++) 
		fprintf(this->velocityFile,"%f ",velocities[n]);

	fprintf(this->velocityFile,"\n");
}











