#include <iostream>
#include "math.h"
#include <fstream>
#include "Molecule.h"
#include "System.h"
#include "lib.h"
#include "omp.h"
const double pi = 3.141592654;
const double boltz = 1.3806e-23;    // Boltzmann's constant (J/K)
double mass = 6.63e-26;     	    // Mass of argon atom (kg)
double diam = 3.66e-10;    	  	    // Effective diameter of argon atom (m)
double density = 2.685e25;		    // Number density of argon at STP (m^-3)

double rand_gauss(long *idum) {
  return sqrt( -2.0*log(1.0 - ran0(idum)) )
	        * cos( 6.283185307 * ran0(idum) );
}

System::System(int N, double T) {
	this->idum = new long[1];
	*this->idum = -1;

	this->N = N;
	this->L = 1e-6;
	this->T = T;

	this->initialize();
	
	printf("System initialized.\n");
}

void System::move() {
	for(int n=0; n< this->N; n++ )
		this->molecules[n]->move(this->tau);
 
	//* Check each particle to see if it strikes a wall
	
	vec xwall = zeros<vec>(2,1);
	vec vw = zeros<vec>(2,1);
	vec direction = zeros<vec>(2,1);
	
	xwall(0) = 0;    xwall(1) = L;   // Positions of walls
	vw(0) = -this->vwall;  vw(1) = this->vwall;  // Velocities of walls
	double stdev = this->mpv/sqrt(2.);
	// Direction of particle leaving wall
	direction(0) = 1;  direction(1) = -1;
	Molecule *molecule;
	this->delta_v.zeros();
	
	#pragma omp parallel for
	for(int n=0; n<this->N; n++) {
		molecule = this->molecules[n];
		//* Test if particle strikes either wall
		int flag = 0;
		if( molecule->r(0) <= 0 )
			flag=1;       // Particle strikes left wall
		else if( molecule->r(0) >= L )
			flag=2;       // Particle strikes right wall

		//* If particle strikes a wall, reset its position
		//  and velocity. Record velocity change.
		if( flag > 0 ) {
			this->wallStrikes(flag-1)++;
			double vyInitial = molecule->v(1);
			//* Reset velocity components as biased Maxwellian,
			//  Exponential dist. in x; Gaussian in y and z
			molecule->v(0) = direction(flag-1)*sqrt(-log(1.-ran0(this->idum))) * this->mpv;
			molecule->v(1) = stdev*rand_gauss(this->idum) + vw(flag-1); // Add wall velocity
			molecule->v(2) = stdev*rand_gauss(this->idum);

			// Time of flight after leaving wall
			double dtr = this->tau*(molecule->r(0)-xwall(flag-1))/(molecule->r(0)-molecule->r_old(0));   
			//* Reset position after leaving wall
			molecule->r(0) = xwall(flag-1) + molecule->v(0)*dtr;
			//* Record velocity change for force measurement
			this->delta_v(flag-1) += (molecule->v(1) - vyInitial);
		}
	}
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell
	int jcell;
	for(int n=0; n<this->numberOfCells; n++ ) {
		col += this->cells[n]->collide();
	}

	return col;
}

void System::step() {
	this->steps += 1;
	this->t += 1;

	//* Move all the particles
	this->move();
	//* Sort the particles into cells
	this->sorter->sort();
	
	this->collisions += this->collide();
}

void System::initialize() {
	this->wallStrikes = zeros<vec>(2,1);
	this->delta_v = zeros<vec>(2,1);

	this->volume = pow(this->L,3);
	this->steps = 0;
	
	this->eff_num = density*this->volume/this->N;
	this->cellsPerDimension = 5;
	this->numberOfCells = this->cellsPerDimension*this->cellsPerDimension*this->cellsPerDimension;

  	printf("Each particle represents %d atoms\n",(int)this->eff_num);
  	this->mfp = this->volume/(sqrt(2.0)*pi*diam*diam*this->N*this->eff_num);
  	printf("System width is %.2f mean free paths\n",L/this->mfp);

  	this->mpv = sqrt(2*boltz*this->T/mass);  // Most probable initial velocity

  	cout << "Enter wall velocity as Mach number: ";
  	double vwall_m;
	cin >> vwall_m;

	this->vwall = vwall_m * sqrt(5./3. * boltz*this->T/mass);
	cout << "Wall velocities are " << -this->vwall << " and " << this->vwall << " m/s" << endl;
	this->initMolecules();
	this->initCells();
	this->collisions = 0;

	// Time step should be so most probable velocity runs through 1/5th of a cell during one timestep
	this->tau = 0.2*(this->L/this->cellsPerDimension)/this->mpv;       // Set timestep tau

	cout << "dt = " << this->tau*1e9 << " ns" << endl;

	this->coeff = 0.5*this->eff_num*pi*diam*diam*this->tau/(this->volume/this->numberOfCells);

	this->sorter = new Sorter(this);
	
	this->t = 0;
}

void System::initCells() {
	this->cells = new Cell*[this->numberOfCells];
	for(int n=0;n<this->numberOfCells;n++) {
		this->cells[n] = new Cell(this);
		this->cells[n]->vr_max = 3*this->mpv;
		this->cells[n]->index = n;
	}
}

void System::initMolecules() {
	this->molecules = new Molecule*[this->N];
	for(int n=0;n<this->N;n++) {
		this->molecules[n] = new Molecule(this);
		this->molecules[n]->mass = this->eff_num;
	}

	this->initPositions();
	this->initVelocities();
}

void System::initPositions() {
	for(int n=0; n<this->N; n++ ) {
		this->molecules[n]->r(0) = L*ran0(this->idum);
		this->molecules[n]->r(1) = L*ran0(this->idum);
		this->molecules[n]->r(2) = L*ran0(this->idum);
	}
}

void System::initVelocities() {
	Molecule *molecule;
	vec v_cm = zeros<vec>(3,1);

	for(int n=0; n<this->N; n++ ) {
		molecule = this->molecules[n];
		
		molecule->v(0) = rand_gauss(this->idum)*sqrt(boltz*this->T/mass);
		molecule->v(1) = rand_gauss(this->idum)*sqrt(boltz*this->T/mass);
		molecule->v(2) = rand_gauss(this->idum)*sqrt(boltz*this->T/mass);

		v_cm += molecule->v;
		molecule->v(1) += this->vwall * 2*(molecule->r(0)/L - 0.5);
  	}

  	v_cm /= this->N;
  	
  	// cout << "vcm=" << v_cm << endl;

  	// Remove any center of mass momentum
  	/*
  	for(int n=0; n<this->N; n++ ) {
  		molecule->v -= v_cm;
  	}
  	*/
}

void System::printPositionsToFile(FILE *file) {
	fprintf(file,"%d\n",this->N);
	fprintf(file,"Random comment that must be here\n");
	
	for(int n=0;n<this->N;n++) {
		fprintf(file,"H %f %f %f\n",1000000*this->molecules[n]->r(0),1000000*this->molecules[n]->r(1),1000000*this->molecules[n]->r(2));
	}
	
}