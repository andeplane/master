#include <system.h>
#include <Image.h>
#include <defines.h>

void System::initialize(CIniFile &ini) {
    read_ini_file(ini);
    Image img;

    char *world_base = "../worlds/";
    char world_file[100];
    char initial_world_file[100];

    strcpy(world_file,world_base);
    strcpy(initial_world_file,world_base);
    strcat(world_file,ini.getstring("world").c_str());
    strcat(initial_world_file,ini.getstring("initial_world").c_str());

    mat world = img.readBMP(world_file);
    mat initial_world = img.readBMP(initial_world_file);

    world_grid = new Grid(world,this);
    initial_world_grid = new Grid(initial_world,this);

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

    init_randoms();
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
    printf("Mean free paths per cell: %.2f \n",min(width/cells_x/mfp,height/cells_y/mfp));
    printf("%d atoms per molecule\n",(int)eff_num);
    printf("%d molecules per cell\n",N/numberOfCells);

    printf("dt = %f\n\n",dt);
}

void System::read_ini_file(CIniFile &ini) {
    N       = ini.getint("N");
    width   = ini.getdouble("width");
    height  = ini.getdouble("height");
    cells_x = ini.getint("cells_x");
    cells_y = ini.getint("cells_y");
    T       = ini.getdouble("T");
    wall_temperature = ini.getdouble("wall_temperature");
    threads = ini.getint("threads");
    acceleration = ini.getdouble("acceleration");
    max_x_acceleration = ini.getdouble("max_x_acceleration");
    density = ini.getdouble("density");
    diam = ini.getdouble("diam");
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

void System::init_randoms() {
    randoms = new Random*[threads];
    for(int i=0;i<threads;i++)
        randoms[i] = new Random(-(i+1));
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
            m->initial_r = m->r;

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
