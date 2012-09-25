#include "StatisticsSampler.h"
#include <armadillo>
#include "Molecule.h"

using namespace arma;

StatisticsSampler::StatisticsSampler(System *system) {
	this->system = system;
	this->temperature = true;
	this->pressure = true;
	this->energy = true;
	this->printVelocities = false;
	this->diffusionConstant = true;

	this->temperatureFile = fopen("temperature.dat","w");
	this->pressureFile    = fopen("pressure.dat","w");
	this->energyFile      = fopen("energy.dat","w");
	this->velocityFile = fopen("velocity.dat","w");
	// this->diffusionFile = fopen("diffusion.dat","w");
}

void StatisticsSampler::sample() {
	this->calculateTemperature();
	this->calculateEnergy();
	this->calculatePressure();
	this->calculateVelocities();
	this->calculateDiffusionConstant();
}

void StatisticsSampler::calculateTemperature() {
	if(!this->temperature) return;
	int N = this->system->N;
	double energy = 0;

	Molecule **molecules = this->system->molecules;
	for(int n=0;n<N;n++) {
		energy += molecules[n]->mass*norm(molecules[n]->v,2);
	}

	energy/=(3*(N-1));

	fprintf(this->temperatureFile, "%f %f \n",this->system->t, energy);
}

void StatisticsSampler::calculateEnergy() {
	if(!this->energy) return;

	int N = this->system->N;
	double E = 0, Ek=0,Ep=0, Ek_temp, Ep_temp;

	Molecule **molecules = this->system->molecules;

	for(int n=0;n<N;n++) {
		Ek_temp = 0.5*molecules[n]->mass*dot(molecules[n]->v,molecules[n]->v);
		Ep_temp = 0;
		Ek += Ek_temp;
		Ep += Ep_temp;
		E += Ek_temp + Ep_temp;
	}

	fprintf(this->energyFile, "%f %f %f %f \n",this->system->t,Ek,Ep,E);
}

void StatisticsSampler::calculatePressure() {
	if(!this->pressure) return;

	// fprintf(this->pressureFile, "%f %f \n",t, this->system->P);
}

void StatisticsSampler::calculateVelocities() {
	if(!this->printVelocities || this->system->steps % 50) return;

	int N = this->system->N;
	
	Molecule **molecules = this->system->molecules;
	
	if(!this->system->steps) fprintf(this->velocityFile, "%d %d %d\n",N,N,N);
	Molecule *molecule;
	double vsq = 0;

	for(int n=0;n<N;n++) {
		molecule = molecules[n];
		
		fprintf(this->velocityFile, "%f %f %f %f \n",this->system->t, molecule->v(0),molecule->v(1),molecule->v(2));
	}
}

void StatisticsSampler::calculateDiffusionConstant() {
	if(!this->diffusionConstant || !this->system->t) return;

/*
	int N = this->system->N;
	Molecule **molecules = this->system->molecules;
	double rsquared = 0;
	for(int n=0;n<N;n++) {
		rsquared += molecules[n]->squaredDistanceFromInitialPosition();
	}
	rsquared /= N;
	double D = rsquared/(6*t);
	
	fprintf(this->diffusionFile, "%f %f\n",t, D);
	// fprintf(this->diffusionConstantFile, "%f %f\n",t, D);
	// fprintf(this->diffusionFile, "Nothing is happening");
	*/
}

