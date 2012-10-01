#include "StatisticsSampler.h"
#include <armadillo>
#include "Molecule.h"
#include "Cell.h"
#include <math.h>
#include "Sorter.h"

using namespace arma;

#define TEMPERATURE 0
#define ENERGY 1
#define PRESSURE 2
#define VELOCITY 3
#define DIFFUSION 4
#define VISCOSITY 5

StatisticsSampler::StatisticsSampler(System *system, int timesteps) {
	this->system = system;
	this->numberOfSamples = new int[6];
	for(int i=0;i<6;i++)
		this->numberOfSamples[i] = 0;
	this->printTemperature = true;

	this->temperatureFile = fopen("temperature.dat","w");

	// Viscosity
	this->delta_v_tot = zeros<vec>(2,1);
	this->delta_v_err = zeros<vec>(2,1);
}

void StatisticsSampler::sample() {
	this->calculateTemperature();
	this->calculateViscosity();
}

void StatisticsSampler::printViscosity() {
	double mass = 6.63e-26;     	    // Mass of argon atom (kg)
	double density = 2.685e25;    // Number density of argon at STP (m^-3)
	double eta = 5.*M_PI/32.*mass*density*(2./sqrt(M_PI)*system->mpv)*system->mfp; // Theoretical value

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
	
	double vgrad = 2*system->vwall/L;    // Velocity gradient

	double visc = 0.5*(-force(0)+force(1))/vgrad;  // Average viscosity
	double viscerr = 0.5*(f_err(0)+f_err(1))/vgrad;  // Error
	cout << "Viscosity = " << visc << " +/- " << viscerr
	   << "N s/m^2" << endl;
	
	cout << "Theoretical value of viscoisty is " << eta
   << "N s/m^2" << endl;
}

void StatisticsSampler::calculateViscosity() {
	System *system = this->system;

	this->delta_v_tot += system->delta_v;
	this->delta_v_err(0) += system->delta_v(0)*system->delta_v(0);
	this->delta_v_err(1) += system->delta_v(1)*system->delta_v(1);

	this->numberOfSamples[VISCOSITY]++;
}

void StatisticsSampler::calculateTemperature() {
	if(!this->printTemperature) return;

	double mass = 6.63e-26;
	double boltz = 1.3806e-23;

	long N = this->system->N;
	double energy = 0;

	Molecule **molecules = this->system->molecules;
	for(int n=0;n<N;n++) {
		energy += molecules[n]->mass*dot(molecules[n]->v,molecules[n]->v);
	}

	// energy/=(3*(N-1));
	double T = energy*mass/(3*N*boltz);
	
	fprintf(this->temperatureFile, "%f %f \n",this->system->t, T);
}
