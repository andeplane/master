#include <iostream>
#include <math.h>
#include <fstream>
#include "Molecule.h"
#include "System.h"
#include <time.h>
#include "defines.h"
#include "CollisionObject.h"

const double boltz = 1.0;    // Boltzmann's constant (J/K)
double mass = 1.0;     	    // Mass of argon atom (kg)
double diam = 3.66e-4;    	  	    // Effective diameter of argon atom (m)
double density = 26850000;		    // Number density of argon at STP (L^-3)


System::System(int N, double T) {
	this->N = N;
	this->L = 1.0;
	this->T = T;
	
	printf("Initializing system...");
	this->initialize();
	printf("done.\n\n");
	printf("%d atoms per molecule\n",(int)this->eff_num);
  	printf("%d cells in each dimension\n",this->cellsPerDimension);
  	printf("%d particles in each cell \n",this->N/this->numberOfCells);
  	printf("%.2f mean free paths (system width)\n",this->L/this->mfp);
  	printf("%.2f mean free paths (cell width)\n",this->L/(this->cellsPerDimension*this->mfp));
  	printf("dt = %f\n\n",this->dt);
}

void System::move() {
	for(int n=0; n< this->N; n++ )
        this->molecules[n]->move(this->dt,randoms[0]);
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell

	int numCells = this->numberOfCells;

	int n;
// #pragma omp parallel
// {
	int local_col = 0;
    // #pragma omp for
	for(n=0; n<numCells; n++ ) {
        local_col += this->cells[n]->collide(randoms[0]);
	}

	col += local_col;
// }


	return col;
}

void System::accelerate() {
    double a = 1;

	for(int n=0;n<this->N;n++) {
		this->molecules[n]->v(0) += a*this->dt;
	}
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

	t0 = clock();

	double E = 0;
	vec p = zeros<vec>(3,1);
	double density = 0;

   //  #pragma omp parallel for
	for(int n=0;n<this->numberOfCells;n++) {
		this->cells[n]->sampleStatistics();
		E += this->cells[n]->energy;
		p += this->cells[n]->momentum;
		density += this->cells[n]->density;
	}

	this->time_consumption[SAMPLE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

	t0 = clock();
	this->collisions += this->collide();
	this->time_consumption[COLLIDE] += ((double)clock()-t0)/CLOCKS_PER_SEC;
	this->accelerate();
}

void System::initialize() {
	this->time_consumption = new double[4];
	for(int i=0;i<4;i++)
		this->time_consumption[i] = 0;

    this->volume = pow(this->L,2);
	this->steps = 0;
    this->eff_num = density*this->volume/this->N;

	this->mfp = this->volume/(sqrt(2.0)*M_PI*diam*diam*this->N*this->eff_num);
	this->mpv = sqrt(this->T);  // Most probable initial velocity
	
	this->cellsPerDimension = 3*this->L/this->mfp;
	this->numberOfCells = this->cellsPerDimension*this->cellsPerDimension;
	
	while(this->N/this->numberOfCells < 20) {
		this->cellsPerDimension--;
		this->numberOfCells = this->cellsPerDimension*this->cellsPerDimension;
	}

	this->dt = 0.2*(this->L/this->cellsPerDimension)/this->mpv;       // Set timestep dt
	this->coeff = 0.5*this->eff_num*M_PI*diam*diam*this->dt/(this->volume/this->numberOfCells);

    this->randoms = new Random*[16];
    for(int i=0;i<16;i++)
        this->randoms[i] = new Random(-(i+1));
	
	this->initMolecules();
	this->initCells();
	this->initWalls();
	this->initObjects();
	this->collisions = 0;

	this->sorter = new Sorter(this);
	
	this->t = 0;
}

void System::initObjects() {
    return;
/*
	this->objects = new CollisionObject*[1];
	vec center = zeros<vec>(2,1);
	center(0) = this->L/2;
	center(1) = this->L/2;

	this->objects[0] = new Box(this, center, this->L/3, this->L/3, this->T);
    */
}

void System::initWalls() {
	this->walls = new Wall*[2];
	this->walls[0] = new Wall(0,this->T,false);
	this->walls[1] = new Wall(L,this->T,true);
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
		this->molecules[n]->atoms = this->eff_num;
	}

	this->initPositions();
	this->initVelocities();
}

void System::initPositions() {

	for(int n=0; n<this->N; n++ ) {
        this->molecules[n]->r(0) = L*randoms[0]->nextDouble();
        this->molecules[n]->r(1) = L*randoms[0]->nextDouble();
	}
}

void System::initVelocities() {
	Molecule *molecule;

	for(int n=0; n<this->N; n++ ) {
		molecule = this->molecules[n];

        molecule->v(0) = randoms[0]->nextGauss()*sqrt(3.0/2*this->T);
        molecule->v(1) = randoms[0]->nextGauss()*sqrt(3.0/2*this->T);
  	}
}

void System::printPositionsToFile(FILE *file) {
	fprintf(file,"%d\n",this->N);
	fprintf(file,"Random comment that must be here\n");

	for(int n=0;n<this->N;n++) {
        fprintf(file,"H %.10f %.10f 0\n",this->molecules[n]->r(0),this->molecules[n]->r(1));
    }
	
}
