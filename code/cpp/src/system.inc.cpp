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

    N    = settings->number_of_particles;
    Lx   = settings->Lx;
    Ly   = settings->Ly;
    Lz   = settings->Lz;
    temperature       = unit_converter->temperature_from_SI(settings->temperature);;
    wall_temperature = unit_converter->temperature_from_SI(settings->wall_temperature);

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
    int c_x,c_y,c_z;

    for(int i=0;i<world_grid->cols;i++) {
        for(int j=0;j<world_grid->rows;j++) {
            for(int k=0;k<world_grid->slices;k++) {
                // Update the effective cell volume. A cell may contain 50% of solid material

                c_x = i/(float)world_grid->cols*settings->cells_x;
                c_y = j/(float)world_grid->rows*settings->cells_y;
                c_z = k/(float)world_grid->slices*settings->cells_z;
                c = cells[c_x][c_y][c_z];
                c->total_pixels++;

                c->pixels += world(j,i) == 0;
                porosity += world(j,i) == 0;
            }
        }
    }


    porosity /= world_grid->rows*world_grid->cols*world_grid->slices;
    volume = Lx*Ly*Lz*porosity;

    eff_num = density*volume/N;
    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*N*eff_num);
    mpv = sqrt(temperature);  // Most probable initial velocity
    double min_cell_size = min( min(Lx/settings->cells_x,Ly/settings->cells_y), Lz/settings->cells_z);

    dt = 0.2*min_cell_size/mpv;       // Set timestep dt
    dt = settings->dt;
    dt *= settings->dt_factor;

    for(int i=0;i<settings->cells_x;i++) {
        for(int j=0;j<settings->cells_y;j++) {
            for(int k=0;k<settings->cells_z;k++) {
                c = cells[i][j][k];
                c->vr_max = 3*mpv;
                c->update_volume();
            }
        }
    }

    init_randoms();

    if(settings->load_previous_state) {
        io->load_state_from_file_binary();
    } else {
        init_molecules();
    }
    int number_of_cells = settings->cells_x*settings->cells_y*settings->cells_z;

    sorter = new Sorter(this);

    printf("done.\n\n");
    printf("%d molecules\n",N);
    printf("%d (%d x %d x %d) cells\n",number_of_cells,settings->cells_x,settings->cells_y,settings->cells_z);
    printf("Porosity: %f\n",porosity);
    printf("System volume: %f\n",volume);
    printf("System size: %.2f x %.2f x %.2f μm\n",Lx,Ly,Lz);
    printf("System size (mfp): %.2f x %.2f x %.2f \n",Lx/mfp,Ly/mfp,Lz/mfp);
    printf("Global Kn: %.2f x %.2f x %.2f \n",mfp/Lx,mfp/Ly, mfp/Lz);
    printf("Mean free paths per cell: %.2f \n",min( min(Lx/settings->cells_x/mfp,Ly/settings->cells_y/mfp), Lz/settings->cells_z/mfp));
    printf("%d atoms per molecule\n",(int)eff_num);
    printf("%d molecules per cell\n",N/number_of_cells);
    printf("Eff_num: %f\n",eff_num);


    printf("dt = %f\n\n",dt);
}

void System::init_cells() {
    cells.resize(settings->cells_x);

    for(int i=0;i<settings->cells_x;i++) {
        vector< vector<Cell*> > &vx = cells[i];
        vx.resize(settings->cells_y);

        for(int j=0;j<settings->cells_y;j++) {
            vector<Cell*> &vy = cells[i][j];
            vy.reserve(settings->cells_z);

            for(int k=0;k<settings->cells_z;k++) {
                Cell *c = new Cell(this);
                c->i = i;
                c->j = j;
                c->k = k;
                c->vr_max = 3*mpv;

                cells[i][j].push_back(c);
            }
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
            m->r[0] = Lx*rnd->nextDouble();
            m->r[1] = Ly*rnd->nextDouble();
            m->r[2] = Lz*rnd->nextDouble();
            did_collide = initial_world_grid->get_grid_point(m->r)->is_wall;
        }

        m->initial_r[0] = m->r[0];
        m->initial_r[1] = m->r[1];
        m->initial_r[2] = m->r[2];
    }
}

void System::init_velocities() {
    Molecule *m;
    for(int n=0; n<N; n++ ) {
        m = molecules[n];
        m->v[0] = rnd->nextGauss()*sqrt(3.0/2*temperature);
        m->v[1] = rnd->nextGauss()*sqrt(3.0/2*temperature);
        m->v[2] = rnd->nextGauss()*sqrt(3.0/2*temperature);
    }
}
