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


System::System(int _N, double _T) {
    N = _N;
    L = 1.0;
    T = _T;
	
	printf("Initializing system...");
    initialize();
	printf("done.\n\n");
    printf("%d atoms per molecule\n",(int)eff_num);
    printf("%d cells in each dimension\n",cellsPerDimension);
    printf("%d particles in each cell \n",N/numberOfCells);
    printf("%.2f mean free paths (system width)\n",L/mfp);
    printf("%.2f mean free paths (cell width)\n",L/(cellsPerDimension*mfp));
    printf("dt = %f\n\n",dt);
}

void System::move() {
    for(int n=0; n< N; n++ )
        molecules[n]->move(dt,randoms[0]);
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell

    int numCells = numberOfCells;

	int n;
// #pragma omp parallel
// {
	int local_col = 0;
    // #pragma omp for
	for(n=0; n<numCells; n++ ) {
        local_col += cells[n]->collide(randoms[0]);
	}

	col += local_col;
// }


	return col;
}

void System::accelerate() {
    double a = 1;

    for(int n=0;n<N;n++) {
        molecules[n]->v(0) += a*dt;
	}
}

void System::step() {
	time_t t0;

    steps += 1;
    t += 1;

	//* Move all the particles
	t0 = clock();
    move();
    time_consumption[MOVE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

	//* Sort the particles into cells
	t0 = clock();
    sorter->sort();
    time_consumption[SORT] += ((double)clock()-t0)/CLOCKS_PER_SEC;

	t0 = clock();

	double E = 0;
	vec p = zeros<vec>(3,1);
	double density = 0;

   //  #pragma omp parallel for
    for(int n=0;n<numberOfCells;n++) {
        cells[n]->sampleStatistics();
        E += cells[n]->energy;
        p += cells[n]->momentum;
        density += cells[n]->density;
	}

    time_consumption[SAMPLE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

	t0 = clock();
    collisions += collide();
    time_consumption[COLLIDE] += ((double)clock()-t0)/CLOCKS_PER_SEC;
    accelerate();
}

void System::initialize() {
    time_consumption = new double[4];
	for(int i=0;i<4;i++)
        time_consumption[i] = 0;

    volume = pow(L,2);
    steps = 0;
    eff_num = density*volume/N;

    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*N*eff_num);
    mpv = sqrt(T);  // Most probable initial velocity
	
    cellsPerDimension = 3*L/mfp;
    numberOfCells = cellsPerDimension*cellsPerDimension;

    while(N/numberOfCells < 20) {
        cellsPerDimension--;
        numberOfCells = cellsPerDimension*cellsPerDimension;
	}

    dt = 0.2*(L/cellsPerDimension)/mpv;       // Set timestep dt
    coeff = 0.5*eff_num*M_PI*diam*diam*dt/(volume/numberOfCells);

    randoms = new Random*[16];
    for(int i=0;i<16;i++)
        randoms[i] = new Random(-(i+1));
	
    initMolecules();
    initCells();
    initWalls();
    initObjects();
    collisions = 0;

    sorter = new Sorter(this);
	
    t = 0;
}

void System::initObjects() {
    return;
/*
    objects = new CollisionObject*[1];
	vec center = zeros<vec>(2,1);
    center(0) = L/2;
    center(1) = L/2;

    objects[0] = new Box(this, center, L/3, L/3, T);
    */
}

void System::initWalls() {
    walls = new Wall*[2];
    walls[0] = new Wall(0,T,false);
    walls[1] = new Wall(L,T,true);
}

void System::initCells() {
    cells = new Cell*[numberOfCells];
    for(int n=0;n<numberOfCells;n++) {
        cells[n] = new Cell(this);
        cells[n]->vr_max = 3*mpv;
        cells[n]->index = n;
	}
}

void System::initMolecules() {
    molecules = new Molecule*[N];
    for(int n=0;n<N;n++) {
        molecules[n] = new Molecule(this);
        molecules[n]->atoms = eff_num;
	}

    initPositions();
    initVelocities();
}

void System::initPositions() {

    for(int n=0; n<N; n++ ) {
        molecules[n]->r(0) = L*randoms[0]->nextDouble();
        molecules[n]->r(1) = L*randoms[0]->nextDouble();
	}
}

void System::initVelocities() {
	Molecule *molecule;

    for(int n=0; n<N; n++ ) {
        molecule = molecules[n];

        molecule->v(0) = randoms[0]->nextGauss()*sqrt(3.0/2*T);
        molecule->v(1) = randoms[0]->nextGauss()*sqrt(3.0/2*T);
  	}
}

void System::printPositionsToFile(FILE *file) {
    fprintf(file,"%d\n",N);
	fprintf(file,"Random comment that must be here\n");

    for(int n=0;n<N;n++) {
        fprintf(file,"H %.10f %.10f 0\n",molecules[n]->r(0),molecules[n]->r(1));
    }
	
}
