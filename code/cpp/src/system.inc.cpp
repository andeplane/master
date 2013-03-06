#include <system.h>
#include <mpi.h>

void System::initialize(Settings *settings_, int myid_) {
    if(myid==0) cout << "Initializing system..." << endl;

    settings = settings_;
    myid = myid_;
    steps = 0;
    collisions = 0;
    t = 0;

    init_randoms();
    io = new DSMC_IO(this);
    unit_converter = new UnitConverter();

    Lx   = settings->Lx;
    Ly   = settings->Ly;
    Lz   = settings->Lz;

    cell_length_x = Lx/(settings->cells_per_node_x*settings->nodes_x);
    cell_length_y = Ly/(settings->cells_per_node_y*settings->nodes_y);
    cell_length_z = Lz/(settings->cells_per_node_z*settings->nodes_z);

    temperature       = unit_converter->temperature_from_SI(settings->temperature);;
    wall_temperature = unit_converter->temperature_from_SI(settings->wall_temperature);

    acceleration = settings->acceleration;
    max_x_acceleration = settings->max_x_acceleration;
    density = settings->density;
    diam = settings->diam;

    if(myid==0) cout << "Loading world..." << endl;
    world_grid = new Grid(settings->ini_file.getstring("world"),this);
    if(myid==0) cout << "Initializing cells..." << endl;
    thread_control.setup(this);

    MPI_Allreduce(&thread_control.porosity,&porosity_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&thread_control.num_particles,&num_particles_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    volume = Lx*Ly*Lz*porosity_global;

    eff_num = density*volume/num_particles_global;
    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*num_particles_global*eff_num);
    mpv = sqrt(temperature);  // Most probable initial velocity

    dt = settings->dt;

    if(myid==0) {
        int number_of_cells = settings->cells_x*settings->cells_y*settings->cells_z;

        printf("done.\n\n");
        printf("%d molecules\n",num_particles_global);
        printf("%d (%d x %d x %d) cells\n",number_of_cells,settings->cells_x,settings->cells_y,settings->cells_z);
        printf("Porosity: %f\n",porosity_global);
        printf("System volume: %f\n",volume);
        printf("System size: %.2f x %.2f x %.2f μm\n",Lx,Ly,Lz);
        printf("System size (mfp): %.2f x %.2f x %.2f \n",Lx/mfp,Ly/mfp,Lz/mfp);
        printf("Global Kn: %.2f x %.2f x %.2f \n",mfp/Lx,mfp/Ly, mfp/Lz);
        printf("Mean free paths per cell: %.2f \n",min( min(Lx/settings->cells_x/mfp,Ly/settings->cells_y/mfp), Lz/settings->cells_z/mfp));
        printf("%d atoms per molecule\n",(int)eff_num);
        printf("%d molecules per cell\n",num_particles_global/number_of_cells);

        printf("dt = %f\n\n",dt);
    }
}

void System::init_randoms() {
    long seed = time(NULL);
    seed = -1;
    rnd = new Random(-seed);
}
