#include <iostream>
#include <math.h>
#include <fstream>
#include <Molecule.h>
#include <System.h>
#include <time.h>
#include <omp.h>
#include <defines.h>

const double boltz = 1.0;    // Boltzmann's constant (J/K)
double mass = 1.0;     	    // Mass of argon atom (kg)
double diam = 3.66e-4;    	  	    // Effective diameter of argon atom (m)
double density = 26850000;		    // Number density of argon at STP (L^-3)

void System::move() {
    Random *rnd = randoms[0];
    for(int n=0; n< N; n++ )
        molecules[n]->move(dt,rnd);
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell

    Random *rnd = randoms[0];
    for(int i=0;i<cells_x;i++)
        for(int j=0;j<cells_y;j++)
        col += cells[i][j]->collide(rnd);

	return col;
}

void System::accelerate() {
    double a = 10;
    for(int n=0;n<N;n++)
        molecules[n]->v(0) += a*dt;
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
    Cell *c;
    for(int i=0;i<cells_x;i++)
        for(int j=0;j<cells_y;j++) {
            c = cells[i][j];
            c->sampleStatistics();
            E += c->energy;
            p += c->momentum;
            density += c->density;
	}

    time_consumption[SAMPLE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

	t0 = clock();
    collisions += collide();
    time_consumption[COLLIDE] += ((double)clock()-t0)/CLOCKS_PER_SEC;
    accelerate();
}

void System::initialize(CIniFile &ini) {
    N       = ini.getint("N");
    width   = ini.getdouble("width");
    height  = ini.getdouble("height");
    cells_x = ini.getint("cells_x");
    cells_y = ini.getint("cells_y");
    T       = ini.getdouble("T");
    threads = ini.getdouble("threads");

    printf("Initializing system...");
    time_consumption = new double[4];
	for(int i=0;i<4;i++)
        time_consumption[i] = 0;

    volume = width*height;
    steps = 0;
    collisions = 0;
    t = 0;
    eff_num = density*volume/N;

    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*N*eff_num);
    mpv = sqrt(T);  // Most probable initial velocity

    numberOfCells = cells_x*cells_y;
    double cell_size = width/cells_x;

    dt = 0.2*cell_size/mpv;       // Set timestep dt
    coeff = 0.5*eff_num*M_PI*diam*diam*dt/(volume/numberOfCells);

    randoms = new Random*[16];
    for(int i=0;i<16;i++)
        randoms[i] = new Random(-(i+1));
	
    initMolecules();
    initCells();
    initWalls();

    sorter = new Sorter(this);

    printf("done.\n\n");
    printf("%d atoms per molecule\n",(int)eff_num);
    printf("%d particles in each cell \n",N/numberOfCells);
    printf("dt = %f\n\n",dt);
}

void System::initWalls() {
    walls = new Wall*[2];
    walls[0] = new Wall(0,T,false);
    walls[1] = new Wall(height,T,true);
}

void System::initCells() {
    cells = new Cell**[cells_x];

    for(int i=0;i<cells_x;i++) {
        cells[i] = new Cell*[cells_y];

        for(int j=0;j<cells_y;j++) {
            cells[i][j] = new Cell(this);
            cells[i][j]->vr_max = 3*mpv;
        }
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
    Molecule *m;
    for(int n=0; n<N; n++ ) {
        m = molecules[n];
        m->r(0) = width*randoms[0]->nextDouble();
        m->r(1) = height*randoms[0]->nextDouble();
	}
}

void System::initVelocities() {
    Molecule *m;
    for(int n=0; n<N; n++ ) {
        m = molecules[n];

        m->v(0) = randoms[0]->nextGauss()*sqrt(3.0/2*T);
        m->v(1) = randoms[0]->nextGauss()*sqrt(3.0/2*T);
  	}
}

void System::printPositionsToFile(FILE *file) {
    fprintf(file,"%d\n",N);
	fprintf(file,"Random comment that must be here\n");

    for(int n=0;n<N;n++) {
        fprintf(file,"H %.10f %.10f 0\n",molecules[n]->r(0),molecules[n]->r(1));
    }
	
}
