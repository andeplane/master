#include <statisticssampler.h>
#include <mdio.h>
#include <settings.h>
#include <unitconverter.h>
#include <system.h>
#include <mpi.h>
#include <mdtimer.h>
#include <iomanip>

StatisticsSampler::StatisticsSampler(System *system_) {
    system = system_;
    settings = system->settings;
    temperature_sampled_at = -1;
    kinetic_energy_sampled_at = -1;
    potential_energy_sampled_at = -1;
    pressure_sampled_at = -1;
}

void StatisticsSampler::sample_momentum_cm() {
    double v_cm_local[3];

    v_cm[0] = 0; v_cm[1] = 0; v_cm[2] = 0;
    v_cm_local[0] = 0; v_cm_local[1] = 0; v_cm_local[2] = 0;

    for(unsigned int i=system->num_atoms_frozen;i<system->num_atoms_local;i++) {
        v_cm_local[0] += system->velocities[3*i+0];
        v_cm_local[1] += system->velocities[3*i+1];
        v_cm_local[2] += system->velocities[3*i+2];
    }

    MPI_Allreduce(v_cm_local,v_cm,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    p_cm[0] = v_cm[0]*settings->mass;
    p_cm[1] = v_cm[1]*settings->mass;
    p_cm[2] = v_cm[2]*settings->mass;
}

void StatisticsSampler::sample_kinetic_energy() {
    if(system->steps == kinetic_energy_sampled_at) return;

    kinetic_energy = 0;
    double kinetic_energy_global = 0;

    for(unsigned int i=system->num_atoms_frozen;i<system->num_atoms_local;i++) {
        kinetic_energy += 0.5*settings->mass*(system->velocities[3*i+0]*system->velocities[3*i+0] + system->velocities[3*i+1]*system->velocities[3*i+1] + system->velocities[3*i+2]*system->velocities[3*i+2]);
    }
    MPI_Allreduce(&kinetic_energy, &kinetic_energy_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    kinetic_energy = kinetic_energy_global;
    kinetic_energy_sampled_at = system->steps;
}

void StatisticsSampler::sample_potential_energy() {
    if(system->steps == potential_energy_sampled_at) return;

    potential_energy = 0;
    MPI_Reduce(&system->potential_energy, &potential_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    potential_energy_sampled_at = system->steps;
}

void StatisticsSampler::sample_temperature() {
    if(system->steps == temperature_sampled_at) return;
    sample_kinetic_energy();
    double kinetic_energy_per_atom = kinetic_energy / system->num_atoms_free_global;
    temperature = 2.0/3*kinetic_energy_per_atom;

    temperature_sampled_at = system->steps;
}

void StatisticsSampler::sample_pressure() {
    if(system->steps == pressure_sampled_at) return;
    sample_temperature();

    pressure = 0;
    MPI_Reduce(&system->pressure_forces,&pressure,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if(system->myid == 0) {
        pressure /= 3*system->volume;
        pressure += system->num_atoms_free_global/system->volume*temperature;
    }

    pressure_sampled_at = system->steps;
}

void StatisticsSampler::sample() {
    system->mdtimer->start_sampling();
    double t_in_pico_seconds = system->unit_converter->time_to_SI(system->t)*1e12;
    sample_temperature();
    sample_potential_energy();
    sample_pressure();
    sample_velocity_distribution();

    if(system->myid == 0) {
        double potential_energy_per_atom = potential_energy/system->num_atoms_free_global;
        double kinetic_energy_per_atom = kinetic_energy/system->num_atoms_free_global;

        fprintf(system->mdio->energy_file, "%.15f %.15f %.15f %.15f %.15f\n",t_in_pico_seconds,
                system->unit_converter->energy_to_ev(kinetic_energy_per_atom),
                system->unit_converter->energy_to_ev(potential_energy_per_atom),
                system->unit_converter->energy_to_ev(kinetic_energy_per_atom+potential_energy_per_atom),
                system->unit_converter->temperature_to_SI(temperature)
                );

        fprintf(system->mdio->pressure_file, "%.15f %.15f\n",t_in_pico_seconds,
                system->unit_converter->pressure_to_SI(pressure)
                );
        cout.setf(ios::fixed);
        cout.precision(5);
        cout << "Timestep " << setw(6) << system->steps << "   t=" << t_in_pico_seconds << " ps   T=" << system->unit_converter->temperature_to_SI(temperature) << " K" << endl;
    }

    system->mdtimer->end_sampling();
}

void StatisticsSampler::sample_velocity_distribution() {
//    int bins_per_dimension = 50;
//    int bins_per_dimension_squared = bins_per_dimension*bins_per_dimension;
//    double *vel = new double[3*bins*bins*bins];

//    for(unsigned int i=system->num_atoms_frozen;i<system->num_atoms_local;i++) {
//        int bin_x = (system->positions[i][0]+system->origo[0]) / system->system_length[0]*bins_per_dimension;
//        int bin_y = (system->positions[i][1]+system->origo[1]) / system->system_length[1]*bins_per_dimension;
//        int bin_z = (system->positions[i][2]+system->origo[2]) / system->system_length[2]*bins_per_dimension;
//        int index = bin_x*bins_per_dimension_squared + bin_y*bins_per_dimension + bin_z;
//        vel[index][0] += system->velocities[3*i+0];
//        vel[index][1] += system->velocities[3*i+1];
//        vel[index][2] += system->velocities[3*i+2];
//    }
//    MPI_Allreduce(&kinetic_energy, &kinetic_energy_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int num_bins_per_dimension = 50;
    int num_bins = num_bins_per_dimension*num_bins_per_dimension;
    double *vel = new double[3*num_bins];
    double *vel_global = new double[3*num_bins];

    int *count = new int[num_bins];
    int *count_global = new int[num_bins];
    memset((void*)count,0,num_bins*sizeof(int));
    memset((void*)count_global,0,num_bins*sizeof(int));
    memset((void*)vel,0,3*num_bins*sizeof(double));
    memset((void*)vel_global,0,3*num_bins*sizeof(double));

    for(unsigned int i=system->num_atoms_frozen;i<system->num_atoms_local;i++) {
        int bin_x = (system->positions[i][0]+system->origo[0]) / system->system_length[0]*num_bins_per_dimension;
        int bin_y = (system->positions[i][1]+system->origo[1]) / system->system_length[1]*num_bins_per_dimension;
        int bin_z = (system->positions[i][2]+system->origo[2]) / system->system_length[2]*num_bins_per_dimension;
        int index = bin_y*num_bins_per_dimension + bin_z;

        vel[3*index+0] += system->velocities[3*i+0];
        vel[3*index+1] += system->velocities[3*i+1];
        vel[3*index+2] += system->velocities[3*i+2];
        count[index]++;
    }

    MPI_Reduce(vel,vel_global,3*num_bins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(count,count_global,num_bins,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(system->myid==0) {
        for(int i=0;i<num_bins;i++) {
            if(count_global[i]>0) {
                vel_global[3*i+0] /= count_global[i];
                vel_global[3*i+1] /= count_global[i];
                vel_global[3*i+2] /= count_global[i];
            }

            fprintf(system->mdio->velocity_file,"%f %f %f ",vel_global[3*i+0],vel_global[3*i+1],vel_global[3*i+2]);
        }
        fprintf(system->mdio->velocity_file,"\n");
    }

    delete vel;
    delete vel_global;
    delete count;
    delete count_global;
}
