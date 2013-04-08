#include <statisticssampler.h>

#include <armadillo>
#include <molecule.h>
#include <cell.h>
#include <math.h>
#include <unitconverter.h>
#include <dsmc_io.h>
#include <system.h>
#include <mpi.h>
#include <threadcontrol.h>

using namespace arma;

StatisticsSampler::StatisticsSampler(System *system_) {
    system = system_;
    settings = system->settings;
    temperature_sampled_at = -1;
    kinetic_energy_sampled_at = -1;
    velocity_distribution_sampled_at = -1;
}

void StatisticsSampler::sample_kinetic_energy() {
    if(system->steps == kinetic_energy_sampled_at) return;
    kinetic_energy_sampled_at = system->steps;

    kinetic_energy = 0;
    double kinetic_energy_global = 0;

    for(unsigned int i=0;i<system->thread_control.num_molecules;i++) {
        kinetic_energy += 0.5*settings->mass*(system->thread_control.v[3*i+0]*system->thread_control.v[3*i+0] + system->thread_control.v[3*i+1]*system->thread_control.v[3*i+1] + system->thread_control.v[3*i+2]*system->thread_control.v[3*i+2]);
    }
    MPI_Allreduce(&kinetic_energy, &kinetic_energy_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    kinetic_energy = kinetic_energy_global;
}

void StatisticsSampler::sample_temperature() {
    if(system->steps == temperature_sampled_at) return;
    temperature_sampled_at = system->steps;

    sample_kinetic_energy();
    double kinetic_energy_per_atom = kinetic_energy / system->num_molecules_global;
    temperature = 2.0/3*kinetic_energy_per_atom;
}

void StatisticsSampler::sample() {
    if(settings->statistics_interval && system->steps % settings->statistics_interval != 0) return;

    double t_in_nano_seconds = system->unit_converter->time_to_SI(system->t)*1e9;

    sample_velocity_distribution();
    sample_temperature();

    long collisions_global = 0;

    MPI_Reduce(&system->collisions, &collisions_global, 1, MPI_LONG,
                   MPI_SUM, 0, MPI_COMM_WORLD);

    if(system->myid == 0) {
        fprintf(system->io->energy_file, "%f %f %f\n",t_in_nano_seconds, system->unit_converter->energy_to_eV(kinetic_energy), system->unit_converter->temperature_to_SI(temperature));
        cout << system->steps << "   t=" << t_in_nano_seconds << "   T=" << system->unit_converter->temperature_to_SI(temperature) << "   Collisions: " <<  collisions_global <<  endl;
    }
}

void StatisticsSampler::sample_velocity_distribution() {
    if(system->steps == velocity_distribution_sampled_at) return;
    velocity_distribution_sampled_at = system->steps;

    int num_bins_per_dimension = 50;
    int num_bins = num_bins_per_dimension*num_bins_per_dimension;
    double *vel = new double[3*num_bins];
    int *count = new int[num_bins];
    memset((void*)count,0,num_bins*sizeof(int));
    memset((void*)vel,0,3*num_bins*sizeof(double));

    for(unsigned int i=0;i<system->thread_control.num_molecules;i++) {
        int bin_x = system->thread_control.r[3*i+0] / system->length[0]*num_bins_per_dimension;
        int bin_y = system->thread_control.r[3*i+1] / system->length[1]*num_bins_per_dimension;
        // int bin_z = system->thread_control.r[3*i+2] / system->length[2]*num_bins_per_dimension;
        int index = bin_x*num_bins_per_dimension + bin_y;
        vel[3*index+0] += system->thread_control.v[3*i+0];
        vel[3*index+1] += system->thread_control.v[3*i+1];
        vel[3*index+2] += system->thread_control.v[3*i+2];
        count[index]++;
    }

    for(int i=0;i<num_bins;i++) {
        if(count[i]>0) {
            vel[3*i+0] /= count[i];
            vel[3*i+1] /= count[i];
            vel[3*i+2] /= count[i];
        }

        fprintf(system->io->velocity_file,"%f %f %f ",vel[3*i+0],vel[3*i+1],vel[3*i+2]);
    }
    fprintf(system->io->velocity_file,"\n");

    delete vel;
    delete count;
}

/*
void StatisticsSampler::calculate_velocity_profile() {
    if(!print_velocity_profile) return;
    if(++velocity_profile_samples % sample_every_n) return;

    int N = 100;
	Molecule *molecule;
    vec velocities = zeros<vec>(N,1);
    vec velocity_count= zeros<vec>(N,1);

	int n;


    for(int i=0;i<system->N;i++) {
        molecule = system->molecules[i];
        n = N*molecule->r(1)/system->height;
        velocities(n) += molecule->v(0);
        velocity_count(n)++;
	}

    for(n=0;n<N;n++) {
        if(velocity_count(n)>0)
            velocities(n) /= velocity_count(n);
    }
	
	for(int n=0;n<N;n++) 
        fprintf(velocity_file,"%f ",velocities(n));

    fprintf(velocity_file,"\n");
}

double StatisticsSampler::get_temperature() {
    if(!print_temperature) return 0;

    return temperature_sum/(temperature_samples/ini->getint("sample_every_n"));
}

void StatisticsSampler::calculate_velocity_field() {
    if(!print_velocity_field) return;

    Cell *c;
    for(int j=0;j<system->settings->cells_y;j++) {
        for(int i=0;i<system->settings->cells_x;i++) {
            c = system->cells[i][j];

            fprintf(velocity_field_file_x,"%.5f ",c->momentum(0));
            fprintf(velocity_field_file_y,"%.5f ",c->momentum(1));
        }

        fprintf(velocity_field_file_x,"\n");
        fprintf(velocity_field_file_y,"\n");
    }
}

double StatisticsSampler::calculate_diffusion_constant() {
    Molecule *molecule;
    double r_squared = 0;
    for(int i=0;i<system->N;i++) {
        molecule = system->molecules[i];
        r_squared += molecule->squared_distance_from_initial_position();
    }
    r_squared /= system->N;

    double diffusion_constant = r_squared/6/system->t;
    return diffusion_constant;
}


mat StatisticsSampler::calculate_global_pressure_tensor() {
    vector<Molecule*> molecules;
    for(int n=0;n<system->N;n++)
        molecules.push_back(system->molecules[n]);

    return calculate_pressure_tensor(molecules);
}

vector<mat> StatisticsSampler::calculate_local_pressure_tensor() {
    vector<mat> pressure_tensors;

    Cell *c;
    for(int i=0;i<system->settings->cells_x;i++) {
        for(int j=0;j<system->settings->cells_y;j++) {
            c = system->cells[i][j];
            vector<Molecule*> molecules;
            for(int n=0;n<c->particles;n++)
                molecules.push_back(system->molecules[c->particle_indices[n]]);
            pressure_tensors.push_back(calculate_pressure_tensor(molecules));
        }
    }

    return pressure_tensors;
}

mat StatisticsSampler::calculate_pressure_tensor(vector<Molecule*> molecules) {
    mat p = zeros<mat>(3,3);
    vec v = zeros<vec>(3,1);
    double rho = 0;
    int N = molecules.size();
    if(!N) return p;

    Molecule *m;

    for(int n=0;n<molecules.size();n++) {
        m = molecules[n];
        v += m->mass*m->v/N;
        rho += m->mass/N;
        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                p(i,j) += m->v(i)*m->v(j)*m->mass/N;
            }
        }
    }
    v /= rho;

    for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            p(i,j) -= v(i)*v(j)*rho;
        }
    }

    return p;
}

void StatisticsSampler::update_cell_statistics() {
    Cell *c; Molecule *m;
    for(int i=0;i<system->settings->cells_x;i++) {
        for(int j=0;j<system->settings->cells_y;j++) {
            c = system->cells[i][j];
            vector<Molecule*> molecules;
            vec momentum = zeros<vec>(3,1);
            double energy = 0;

            for(int n=0;n<c->particles;n++) {
                m = system->molecules[c->particle_indices[n]];

                energy += 0.5*dot(m->v,m->v)*m->mass*m->atoms;
                momentum += m->v*m->mass*m->atoms;

                molecules.push_back(m);
            }
            c->update_pressure_tensor(calculate_pressure_tensor(molecules));
            c->update_momentum(momentum);
            c->update_energy(energy);
            c->update_temperature(2*energy/(3*c->particles));
        }
    }
}
*/
