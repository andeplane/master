#include <system.h>
#include <Image.h>
#include <defines.h>
#include <dsmc_io.h>

void System::initialize(Settings *settings_) {
    settings = settings_;

    printf("Initializing system...");
    io = new DSMC_IO(this);
    unit_converter = new UnitConverter();

    time_consumption = new double[4];
    Image img;
    char *world_base = "../worlds/";
    char world_file[100];
    char initial_world_file[100];

    N       = settings->number_of_particles;
    width   = settings->width;
    height  = settings->height;

    temperature       = unit_converter->temperature_from_SI(settings->temperature);;
    wall_temperature = unit_converter->temperature_from_SI(settings->wall_temperature);
    cout << "Temperature: " << temperature << endl;
    cout << "Temperature wall: " << wall_temperature << endl;

    acceleration = settings->acceleration;
    max_x_acceleration = settings->max_x_acceleration;
    density = settings->density;
    diam = settings->diam;

    steps = 0;
    collisions = 0;
    t = 0;

    for(int i=0;i<4;i++) {
        time_consumption[i] = 0;
    }

    strcpy(world_file,world_base);
    strcpy(initial_world_file,world_base);
    strcat(world_file,settings->ini_file.getstring("world").c_str());
    strcat(initial_world_file,settings->ini_file.getstring("initial_world").c_str());

    mat world = img.readBMP(world_file);
    mat initial_world = img.readBMP(initial_world_file);

    world_grid = new Grid(world,this);
    initial_world_grid = new Grid(initial_world,this);

    init_cells();

    double porosity = 0;
    Cell *c;
    int c_x,c_y;
    for(int i=0;i<world.n_cols;i++) {
        for(int j=0;j<world.n_rows;j++) {
            // Update the effective cell volume. A cell may contain 50% of solid material

            c_x = i/(float)world.n_cols*settings->cells_x;
            c_y = j/(float)world.n_rows*settings->cells_y;
            c = cells[c_x][c_y];
            c->total_pixels++;

            c->pixels += world(j,i) == 0;
            porosity += world(j,i) == 0;
        }
    }

    for(int i=0;i<settings->cells_x;i++) {
        for(int j=0;j<settings->cells_y;j++) {
            c = cells[i][j];
            c->update_volume();
        }
    }

    porosity /= world.n_cols*world.n_rows;
    volume = width*height*porosity;

    eff_num = density*volume/N;

    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*N*eff_num);
    mpv = sqrt(temperature);  // Most probable initial velocity

    double cell_size = width/settings->cells_x;

    dt = 0.2*cell_size/mpv;       // Set timestep dt

    dt *= settings->dt_factor;
    int number_of_cells = settings->cells_x*settings->cells_y;

    coeff = 0.5*eff_num*M_PI*diam*diam*dt/(volume/number_of_cells);

    init_randoms();

    if(settings->load_previous_state) {
        io->load_state_from_file_binary();
    } else {
        init_molecules();
    }

    sorter = new Sorter(this);

    printf("done.\n\n");
    printf("%d molecules\n",N);
    printf("%d (%d x %d) cells\n",number_of_cells,settings->cells_x,settings->cells_y);
    printf("Porosity: %f\n",porosity);
    printf("System volume: %f\n",volume);
    printf("System size: %.2f x %.2f \n",width,height);
    printf("System size (mfp): %.2f x %.2f \n",width/mfp,height/mfp);
    printf("Global Kn: %.2f x %.2f \n",mfp/width,mfp/height);
    printf("Mean free paths per cell: %.2f \n",min(width/settings->cells_x/mfp,height/settings->cells_y/mfp));
    printf("%d atoms per molecule\n",(int)eff_num);
    printf("%d molecules per cell\n",N/number_of_cells);

    printf("dt = %f\n\n",dt);
}

void System::init_cells() {
    cells.resize(settings->cells_x);

    for(int i=0;i<settings->cells_x;i++) {
        ((vector<Cell*>)cells[i]).reserve(settings->cells_y);

        for(int j=0;j<settings->cells_y;j++) {
            Cell *c = new Cell(this);
            c->i = i;
            c->j = j;
            c->vr_max = 3*mpv;

            cells[i].push_back(c);
        }
    }
}

void System::init_randoms() {
    long seed = time(NULL);
    rnd = new Random(-seed);
}

void System::init_molecules() {
    molecules.reserve(N);
    positions = new double[3*N];
    velocities = new double[3*N];
    initial_positions = new double[3*N];
    Molecule *m;
    for(int n=0;n<N;n++) {
        m = new Molecule(this);
        m->atoms = eff_num;
        m->index = n;
        m->r = &positions[3*n];
        m->v = &velocities[3*n];
        m->initial_r = &initial_positions[3*n];
        molecules.push_back(m);
    }

    init_positions();
    init_velocities();
}

void System::init_positions() {
    Molecule *m;

    bool did_collide;

    for(int n=0; n<N; n++ ) {
        did_collide = true;
        m = molecules[n];
        while(did_collide) {
            m->r[0] = width*rnd->nextDouble();
            m->r[1] = width*rnd->nextDouble();
            m->r[2] = 0;
            m->initial_r[0] = m->r[0];
            m->initial_r[1] = m->r[1];
            m->initial_r[2] = m->r[2];

            did_collide = initial_world_grid->get_grid_point(m->r)->is_wall;
        }
    }
}

void System::init_velocities() {
    Molecule *m;
    for(int n=0; n<N; n++ ) {
        m = molecules[n];
        m->v[0] = rnd->nextGauss()*sqrt(2.0/2*temperature);
        m->v[1] = rnd->nextGauss()*sqrt(2.0/2*temperature);
        m->v[2] = 0;
    }
}
