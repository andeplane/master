#include <iostream>
#include <math.h>
#include <fstream>
#include <molecule.h>
#include <system.h>
#include <time.h>
#include <omp.h>
#include <defines.h>

const double boltz = 1.0;    // Boltzmann's constant (J/K)
double mass = 1.0;     	    // Mass of argon atom (kg)
double diam = 4.17e-4;    	  	    // Effective diameter of argon atom (m)
double density = 26850000;		    // Number density of argon at STP (L^-3)
mat readBMP(char* filename);

void System::initialize(CIniFile &ini) {
    N       = ini.getint("N");
    width   = ini.getdouble("width");
    height  = ini.getdouble("height");
    cells_x = ini.getint("cells_x");
    cells_y = ini.getint("cells_y");
    T       = ini.getdouble("T");
    wall_temperature = ini.getdouble("wall_temperature");
    threads = ini.getint("threads");
    mat world   = readBMP((char*)ini.getstring("world").c_str());
    mat initial_world = readBMP((char*)ini.getstring("initial_world").c_str());

    world_grid = new Grid(world,this);
    initial_world_grid = new Grid(initial_world,this);

    acceleration = ini.getdouble("acceleration");
    max_x_acceleration = ini.getdouble("max_x_acceleration");

    printf("Initializing system...");
    time_consumption = new double[4];
    for(int i=0;i<4;i++)
        time_consumption[i] = 0;
    double porosity = 0;
    for(int i=0;i<world.n_cols;i++)
        for(int j=0;j<world.n_rows;j++) {
            porosity += world(j,i) == 0;
        }
    porosity /= world.n_cols*world.n_rows;

    volume = width*height*porosity;

    steps = 0;
    collisions = 0;
    t = 0;
    eff_num = density*volume/N;

    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*N*eff_num);
    mpv = sqrt(T);  // Most probable initial velocity

    numberOfCells = cells_x*cells_y;
    double cell_size = width/cells_x;

    dt = 0.2*cell_size/mpv;       // Set timestep dt
    dt *= ini.getdouble("dt_factor");

    coeff = 0.5*eff_num*M_PI*diam*diam*dt/(volume/numberOfCells);

    randoms = new Random*[16];
    for(int i=0;i<16;i++)
        randoms[i] = new Random(-(i+1));
    initMolecules();
    initCells();
    sorter = new Sorter(this);

    printf("done.\n\n");
    printf("%d molecules\n",N);
    printf("%d (%d x %d) cells\n",numberOfCells,cells_x,cells_y);
    printf("Porosity: %f\n",porosity);
    printf("System volume: %f\n",volume);
    printf("System size: %.2f x %.2f \n",width,height);
    printf("System size (mfp): %.2f x %.2f \n",width/mfp,height/mfp);
    printf("Global Kn: %.2f x %.2f \n",mfp/width,mfp/height);
    printf("Mean free paths per cell: %.2f \n",min(mfp/(width/cells_x),mfp/(height/cells_y)));
    printf("%d atoms per molecule\n",(int)eff_num);
    printf("%d molecules per cell\n",N/numberOfCells);

    printf("dt = %f\n\n",dt);
}

void System::step() {
    time_t t0;

    steps += 1;
    t += 1;

    accelerate();

    //* Move all the particles
    t0 = clock();
    move();
    time_consumption[MOVE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

    //* Sort the particles into cells
    t0 = clock();
    // Arrange all particles in correct cells
    sorter->sort_system();

    time_consumption[SORT] += ((double)clock()-t0)/CLOCKS_PER_SEC;

    t0 = clock();

    collisions += collide();
    time_consumption[COLLIDE] += ((double)clock()-t0)/CLOCKS_PER_SEC;

    t0 = clock();

    double E = 0;
    vec p(3,1);
    double density = 0;

    Cell *c;
    #pragma omp parallel for private(c) num_threads(threads)
    for(int i=0;i<cells_x;i++)
        for(int j=0;j<cells_y;j++) {
            c = cells[i][j];
            c->sampleStatistics();
            // E += c->energy;
            // p += c->momentum;

            // density += c->density;
    }

    time_consumption[SAMPLE] += ((double)clock()-t0)/CLOCKS_PER_SEC;
}

void System::move() {
#pragma omp parallel num_threads(threads)
    {
    Random *rnd = randoms[omp_get_thread_num()];
    #pragma omp for
    for(int n=0; n<N; n++ ) {
        molecules[n]->move(dt,rnd);
    }
    }
}

int System::collide() {
	int col = 0;          // Count number of collisions
	
	//* Loop over cells and process collisions in each cell
    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        Random *rnd = randoms[thread_id];
        int local_col = 0;
        Cell **local_cells = load_balanced_cell_list[thread_id];

        #pragma omp for
        for(int i=0;i<cells_x;i++)
            for(int j=0;j<cells_y;j++)
                local_col += cells[i][j]->collide(rnd);

        #pragma omp critical
        {
           col += local_col;
        }
    }

	return col;
}

void System::accelerate() {
    for(int n=0;n<N;n++)
        if(molecules[n]->r(0) < max_x_acceleration)
            molecules[n]->v(0) += acceleration*dt;
}

void System::initCells() {
    cells = new Cell**[cells_x];

    load_balanced_cell_list = new Cell**[threads];
    cells_in_list = new int[threads];
    for(int i=0;i<threads;i++) {
        // Create cell lists that theoretically can contain all cells
        load_balanced_cell_list[i] = new Cell*[cells_x*cells_y];
        cells_in_list[i] = 0;
    }

    for(int i=0;i<cells_x;i++) {
        cells[i] = new Cell*[cells_y];

        for(int j=0;j<cells_y;j++) {
            cells[i][j] = new Cell(this);
            cells[i][j]->i = i;
            cells[i][j]->j = j;
            cells[i][j]->vr_max = 3*mpv;
        }
	}
}

void System::initMolecules() {
    molecules = new Molecule*[N];
    for(int n=0;n<N;n++) {
        molecules[n] = new Molecule(this);
        molecules[n]->atoms = eff_num;
        molecules[n]->index = n;
        if(n==0)
            molecules[n]->information_carrier = 1;
    }

    initPositions();
    initVelocities();
}

void System::initPositions() {
    Molecule *m;

    bool didCollide;

    for(int n=0; n<N; n++ ) {
        didCollide = true;
        m = molecules[n];
        while(didCollide) {
            m->r(0) = width*randoms[0]->nextDouble();
            m->r(1) = height*randoms[0]->nextDouble();

            didCollide = initial_world_grid->get_grid_point(m->r)->is_wall;
        }
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
        // fprintf(file,"%s %.10f %.10f 0\n",molecules[n]->type,molecules[n]->r(0),molecules[n]->r(1));
        fprintf(file,"%s %.10f %.10f 0\n",molecules[n]->information_carrier ? "O" : "H",molecules[n]->r(0),molecules[n]->r(1));
    }
	
}

mat readBMP(char* filename)
{
    int i;
    FILE* f = fopen(filename, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int fileSize = *(int*)&info[2];

    int width  = *(int*)&info[18];
    int height = *(int*)&info[22];

    int pixels = width*height;

    int bytesLeft = fileSize-54; // 54 bytes for the header, we need the rest

    unsigned char* data = new unsigned char[bytesLeft];
    fread(data, sizeof(unsigned char), bytesLeft, f); // read the rest of the data at once
    fclose(f);

    mat img = zeros<mat>(width,height);

    int pixelCount = 0;
    int pixelIndex = 0;

    int row = 0; // BMP-files start with the lower left pixel, row for row
    int col = 0;

    for(i = 0; i < pixels; i++)
    {
        int r = data[pixelIndex+2]; // RGB values are interchanged in the file format
        int g = data[pixelIndex+1];
        int b = data[pixelIndex];

        double avg = 1.0*(r+g+b)/3.0/255.0; // If we have black/white only, the average is 0 (black) to 255 (white)
        img(col,row) = avg;

        pixelCount++;
        pixelIndex+=3; // Each pixel has 3 bytes, RGB

        if(++col == width) {
            // BMP-format is stupid since it wants each row to have 4n bytes, so it adds
            // the remaining pixels before next row
            int padding = col % 4;

            col = 0;
            pixelIndex+=padding;
            row++;
        }
    }

    return img.t();
}
