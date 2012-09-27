#include "StatisticsSampler.h"
#include <armadillo>
#include "Molecule.h"
#include "Cell.h"

using namespace arma;

#define TEMPERATURE 0
#define ENERGY 1
#define PRESSURE 2
#define VELOCITY 3
#define DIFFUSION 4
#define VISCOSITY 5

StatisticsSampler::StatisticsSampler(System *system) {
	this->system = system;
	this->numberOfSamples = new int[6];
	for(int i=0;i<6;i++)
		this->numberOfSamples[i] = 0;
	this->temperature = true;
	this->pressure = false;
	this->energy = false;
	this->printVelocities = false;
	this->diffusionConstant = false;

	this->temperatureFile = fopen("temperature.dat","w");
	this->pressureFile    = fopen("pressure.dat","w");
	this->energyFile      = fopen("energy.dat","w");
	this->velocityFile = fopen("velocity.dat","w");
	// this->diffusionFile = fopen("diffusion.dat","w");

	this->delta_v_tot = zeros<vec>(2,1);
	this->delta_v_err = zeros<vec>(2,1);
}

void StatisticsSampler::sample() {
	this->calculateTemperature();
	this->calculateEnergy();
	this->calculatePressure();
	this->calculateVelocities();
	this->calculateDiffusionConstant();
}

void StatisticsSampler::printViscosity() {
	double mass = 6.63e-26;     	    // Mass of argon atom (kg)
	double pi = 3.141592654;
	double density = 2.685e25;    // Number density of argon at STP (m^-3)
	double eta = 5.*pi/32.*mass*density*(2./sqrt(pi)*system->mpv)*system->mfp; // Theoretical value

	double L = system->L;
	double t = this->numberOfSamples[VISCOSITY]*system->tau;
	
	vec force = zeros<vec>(2,1);
	vec f_err = zeros<vec>(2,1);

	int n_minus_1 = this->numberOfSamples[VISCOSITY]-1;

	this->delta_v_err(0) = this->delta_v_err(0)/n_minus_1 - this->delta_v_tot(0)/n_minus_1*this->delta_v_tot(0)/n_minus_1;
	this->delta_v_err(1) = this->delta_v_err(1)/n_minus_1 - this->delta_v_tot(1)/n_minus_1*this->delta_v_tot(1)/n_minus_1;
	this->delta_v_err = sqrt(this->delta_v_err*this->numberOfSamples[VISCOSITY]);

	force = system->eff_num*mass*this->delta_v_tot/(t*L*L);
	f_err = system->eff_num*mass*this->delta_v_err/(t*L*L);
	
	cout << "L=" << L << endl;
	cout << "tau=" << system->tau << endl;
	cout << "t=" << t << endl;
	cout << "force(0)=" << force(0) << endl;
	cout << "force(1)=" << force(1) << endl;
	cout << "wwall=" << system->vwall;

	cout << "Force per unit area is" << endl;
	cout << "Left wall:  " << force(0) << " +/- " << f_err(0) << endl;
	cout << "Right wall: " << force(1) << " +/- " << f_err(1) << endl;
	double vgrad = 2*system->vwall/L;    // Velocity gradient
	double visc = 0.5*(-force(0)+force(1))/vgrad;  // Average viscosity
	double viscerr = 0.5*(f_err(0)+f_err(1))/vgrad;  // Error
	cout << "Viscosity = " << visc/4 << " +/- " << viscerr
	   << "N s/m^2" << endl;
	
	cout << "Theoretical value of viscoisty is " << eta
   << "N s/m^2" << endl;
}

void StatisticsSampler::calculateViscosity() {
	System *system = this->system;


	this->numberOfSamples[VISCOSITY]++;
	this->delta_v_tot += system->delta_v;
	this->delta_v_err(0) += system->delta_v(0)*system->delta_v(0);
	this->delta_v_err(1) += system->delta_v(1)*system->delta_v(1);
}

void StatisticsSampler::calculateTemperature() {
	if(!this->temperature) return;

	double mass = 6.63e-26;
	double boltz = 1.3806e-23;

	int N = this->system->N;
	double energy = 0;

	Molecule **molecules = this->system->molecules;
	for(int n=0;n<N;n++) {
		energy += molecules[n]->mass*norm(molecules[n]->v,2);
	}

	// energy/=(3*(N-1));
	double T = 2*energy*mass/(3*N*N*boltz);

	fprintf(this->temperatureFile, "%f %f \n",this->system->t, T);
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

	Cell **cells = this->system->cells;
	double density = 2.685e25;    // Number density of argon at STP (m^-3)
	double boltz = 1.3806e-23;    // Boltzmann's constant (J/K)

	double P = 0;

	double delta_t;
	double delta_v;
	for(int n=0;n<this->system->numberOfCells;n++) {
		Cell *cell = cells[n];
		delta_t = this->system->t - cell->timeForPressureReset;
		delta_v = cell->delta_v;
		P += delta_v/delta_t;
		cell->resetPressureCalculation();
	}

	P *= this->system->molecules[0]->mass;
	// P /= 3*this->system->volume;
	P += density*boltz*this->system->T;

	// printf("I have pressure %f\n",P);


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

