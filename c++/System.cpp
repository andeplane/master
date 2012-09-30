#include <iostream>
#include <math.h>
#include <fstream>
#include "Molecule.h"
#include "System.h"
#include "lib.h"
#include "omp.h"
#include <time.h>
#include "defines.h"

const double boltz = 1.3806e-23;    // Boltzmann's constant (J/K)
double mass = 6.63e-26;     	    // Mass of argon atom (kg)
double diam = 3.66e-10;    	  	    // Effective diameter of argon atom (m)
double density = 2.685e25;		    // Number density of argon at STP (m^-3)

double rand_gauss(long *idum) {
  return sqrt( -2.0*log(1.0 - ran0(idum)) )
	        * cos( 6.283185307 * ran0(idum) );
}

System::System(int N, double T) {
	this->idums = new long*[4];
	for(int i=0;i<4;i++) {
		this->idums[i] = new long[1];
		*this->idums[i] = -(i+1);
	}
	
	this->idum = this->idums[0];

	this->N = N;
	this->L = 1e-6;
	this->T = T;

	this->initialize();
	double E = 0;

	for(int i=0;i<N;i++) {
		E+= 0.5*this->molecules[i]->mass*dot(this->molecules[i]->v,this->molecules[i]->v);
	}

	E *= mass;
	E /= this->volume/this->numberOfCells;
	
	printf("System initialized.\n");
}

void System::move() {
#pragma omp parallel for
	for(int n=0; n< this->N; n++ )
		this->molecules[n]->move(this->tau);
 
	//* Check each particle to see if it strikes a wall
	
	this->delta_v.zeros();
	
#pragma omp parallel
{
	vec xwall = zeros<vec>(2,1);
	vec vw = zeros<vec>(2,1);
	vec direction = zeros<vec>(2,1);
	vec delta_v = zeros<vec>(2,1);
	vec wall_strikes = zeros<vec>(2,1);

	double tau = this->tau;
	xwall(0) = 0;    xwall(1) = L;   // Positions of walls
	vw(0) = -this->vwall;  vw(1) = this->vwall;  // Velocities of walls
	double stdev = this->mpv/sqrt(2.);
	// Direction of particle leaving wall
	direction(0) = 1;  direction(1) = -1;
	int particleCount = this->N;
	long *idum = this->idums[omp_get_thread_num()];
	Molecule *molecule;
	double vyInitial;

#pragma omp for
	for(int n=0; n<particleCount; n++) {
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
			wall_strikes(flag-1)++;
			vyInitial = molecule->v(1);
			//* Reset velocity components as biased Maxwellian,
			//  Exponential dist. in x; Gaussian in y and z
			molecule->v(0) = direction(flag-1)*sqrt(-log(1.-ran0(idum))) * this->mpv;
			molecule->v(1) = stdev*rand_gauss(idum) + vw(flag-1); // Add wall velocity
			molecule->v(2) = stdev*rand_gauss(idum);

			// Time of flight after leaving wall
			double dtr = tau*(molecule->r(0)-xwall(flag-1))/(molecule->r(0)-molecule->r_old(0));   
			//* Reset position after leaving wall
			molecule->r(0) = xwall(flag-1) + molecule->v(0)*dtr;
			//* Record velocity change for force measurement
			delta_v(flag-1) += (molecule->v(1) - vyInitial);
		}
	}
	// Update 'global' variables
	this->wallStrikes += wall_strikes;
	this->delta_v += delta_v;
}
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell
	int jcell;
	int numCells = this->numberOfCells;
	
	System *system = this;
	int n;
	
	#pragma omp parallel for
	for(n=0; n<numCells; n++ ) {
		this->cells[n]->collide();
	}

	return col;
}

void System::step() {
	time_t t0;

	this->steps += 1;
	this->t += 1;

	//* Move all the particles
	t0 = clock();
	this->move();
	this->time_consumption[MOVE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

	//* Sort the particles into cells
	t0 = clock();
	this->sorter->sort();
	this->time_consumption[SORT] += ((double)clock()-t0)/CLOCKS_PER_SEC;
	int cellsWithZero = 0;

	t0 = clock();

	double E = 0;
	vec p = zeros<vec>(3,1);
	double density = 0;

	#pragma omp parallel for
	for(int n=0;n<this->numberOfCells;n++) {
		this->cells[n]->sampleStatistics();
		E += this->cells[n]->energy;
		p += this->cells[n]->momentum;
		density += this->cells[n]->density;
	}

	this->time_consumption[SAMPLE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

	// cout << "Energy: " << E << endl;
	// cout << "Momentum: " << p << endl;
	// cout << "Density: " << density << endl;
	t0 = clock();
	this->collisions += this->collide();
	this->time_consumption[COLLIDE] += ((double)clock()-t0)/CLOCKS_PER_SEC;
}

void System::initialize() {
	this->time_consumption = new double[4];
	for(int i=0;i<4;i++)
		this->time_consumption[i] = 0;

	this->wallStrikes = zeros<vec>(2,1);
	this->delta_v = zeros<vec>(2,1);

	this->volume = pow(this->L,3);
	this->steps = 0;
	this->eff_num = density*this->volume/this->N;
	
	this->mfp = this->volume/(sqrt(2.0)*M_PI*diam*diam*this->N*this->eff_num);
	this->mpv = sqrt(2*boltz*this->T/mass);  // Most probable initial velocity
	// Time step should be so most probable velocity runs through 1/5th of a cell during one timestep

	
	this->cellsPerDimension = 3*this->L/this->mfp;
	this->numberOfCells = this->cellsPerDimension*this->cellsPerDimension*this->cellsPerDimension;
	
	while(this->N/this->numberOfCells < 20) {
		this->cellsPerDimension--;
		this->numberOfCells = this->cellsPerDimension*this->cellsPerDimension*this->cellsPerDimension;
	}

	// this->cellsPerDimension = 15; // Approx number of mean free paths 
	

	this->tau = 0.2*(this->L/this->cellsPerDimension)/this->mpv;       // Set timestep tau
	this->coeff = 0.5*this->eff_num*M_PI*diam*diam*this->tau/(this->volume/this->numberOfCells);



  	cout << "Enter wall velocity as Mach number: ";
  	double vwall_m;
	cin >> vwall_m;

	this->vwall = vwall_m * sqrt(5./3. * boltz*this->T/mass);
	
	this->initMolecules();
	this->initCells();
	this->collisions = 0;

	printf("Each particle represents %d atoms\n",(int)this->eff_num);
  	printf("System width is %.2f mean free paths\n",this->L/this->mfp);
  	printf("The system consists of %d cells in each dimension\n",this->cellsPerDimension);
  	printf("Each cell contains approx. %d particles \n",this->N/this->numberOfCells);
  	cout << "Wall velocities are " << -this->vwall << " and " << this->vwall << " m/s" << endl;
	cout << "dt = " << this->tau*1e9 << " ns" << endl;

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

	for(int n=0; n<this->N; n++ ) {
		molecule = this->molecules[n];
		
		molecule->v(0) = rand_gauss(this->idum)*sqrt(boltz*this->T/mass);
		molecule->v(1) = rand_gauss(this->idum)*sqrt(boltz*this->T/mass);
		molecule->v(2) = rand_gauss(this->idum)*sqrt(boltz*this->T/mass);

		// molecule->v(1) += this->vwall * 2*(molecule->r(0)/L - 0.5);
  	}
}

void System::printPositionsToFile(FILE *file) {
	fprintf(file,"%d\n",this->N);
	fprintf(file,"Random comment that must be here\n");
	
	for(int n=0;n<this->N;n++) {
		fprintf(file,"H %f %f %f\n",1000000*this->molecules[n]->r(0),1000000*this->molecules[n]->r(1),1000000*this->molecules[n]->r(2));
	}
	
}