#include <system.h>
#include <mpi.h>
#include <dsmctimer.h>
#include <moleculemover.h>
#include <settings.h>

void System::initialize(Settings *settings_, int myid_) {
    myid = myid_;
    timer = new DSMCTimer();
    timer->start_system_initialize();

    if(myid==0) cout << "Initializing system..." << endl;

    settings = settings_;
    steps = 0;
    collisions = 0;
    t = 0;

    init_randoms();
    unit_converter = new UnitConverter();

    length[0] = settings->Lx;
    length[1] = settings->Ly;
    length[2] = settings->Lz;

    reservoir_size = length[settings->gravity_direction]*settings->reservoir_fraction/2;

    cell_length_x = length[0]/(settings->cells_per_node_x*settings->nodes_x);
    cell_length_y = length[1]/(settings->cells_per_node_y*settings->nodes_y);
    cell_length_z = length[2]/(settings->cells_per_node_z*settings->nodes_z);

    cells_x = settings->nodes_x*settings->cells_per_node_x;
    cells_y = settings->nodes_y*settings->cells_per_node_y;
    cells_z = settings->nodes_z*settings->cells_per_node_z;
    num_cells_vector[0] = cells_x;
    num_cells_vector[1] = cells_y;
    num_cells_vector[2] = cells_z;

    temperature      = unit_converter->temperature_from_SI(settings->temperature);;
    wall_temperature = unit_converter->temperature_from_SI(settings->wall_temperature);

    density = settings->density;
    diam = settings->diam;

    io = new DSMC_IO(this);

    if(myid==0) cout << "Loading world..." << endl;
    world_grid = new Grid(settings->ini_file.getstring("world"),this);
    if(myid==0) cout << "Initializing thread system..." << endl;

    thread_control.setup(this);

    MPI_Allreduce(&thread_control.porosity,&porosity_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&thread_control.num_molecules,&num_molecules_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    porosity_global /= thread_control.num_nodes;

    volume = length[0]*length[1]*length[2]*porosity_global;
    eff_num = density*volume/num_molecules_global;

    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*num_molecules_global*eff_num);
    mpv = sqrt(temperature);  // Most probable initial velocity

    dt = settings->dt;

    if(myid==0) cout << "Updating cell volume..." << endl;
    thread_control.update_cell_volume();
    mover = new MoleculeMover();
    mover->initialize(this);

    if(myid==0) {
        int number_of_cells = thread_control.my_cells.size();
        int number_of_cells_all = cells_x*cells_y*cells_z;

        printf("done.\n\n");
        printf("%d molecules\n",num_molecules_global);
        printf("%d (%d inactive) cells\n",number_of_cells,number_of_cells_all - number_of_cells);
        printf("Porosity: %f\n",porosity_global);
        printf("System volume: %f\n",length[0]*length[1]*length[2]);
        printf("Effective system volume: %f\n",volume);
        printf("Mean free path: %.4f \n",mfp);
        printf("Mean free paths per cell: %.2f \n",min( min(length[0]/cells_x/mfp,length[1]/cells_y/mfp), length[2]/cells_z/mfp));
        printf("%ld atoms per molecule\n",(unsigned long)eff_num);
        printf("%d molecules per active cell\n",num_molecules_global/number_of_cells);

        printf("dt = %f\n\n",dt);
    }

    timer->end_system_initialize();
}

void System::init_randoms() {
    long seed = time(NULL);
    seed = 1 + myid;
    rnd = new Random(-seed);
}
