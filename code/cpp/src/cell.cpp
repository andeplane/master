#include "cell.h"
#include <math.h>
#include <time.h>
#include <molecule.h>

using namespace std;

Cell::Cell(System *_system) {
    system = _system;
    vr_max = 0;
    particles = 0;
    average_over = 100;
    momentum_time_steps = 0;
    temperature = system->T;
    momentum = zeros<vec>(3,1);
    momentum_change = zeros<vec>(3,1);
    pressure_tensor = zeros<mat>(3,3);
    energy = 0;
    density = 0;
    pixels = 0;
    total_pixels = 0;

    particle_capacity = 100;
    particle_indices = new unsigned int[particle_capacity];
}

bool Cell::cmp(Cell *c1, Cell *c2) {
    return c1->collision_pairs < c2->collision_pairs;
}

void Cell::resize(int n) {
    delete [] particle_indices;
    particle_capacity = n;
    particle_indices = new unsigned int[n];
}

void Cell::reset() {
    particles = 0;
}

void Cell::update_pressure_tensor(mat updated_pressure_tensor) {
    pressure_tensor = ((average_over-1)/average_over)*pressure_tensor+1.0/average_over*updated_pressure_tensor;
}

void Cell::update_momentum(vec updated_momentum) {

    momentum = ((average_over-1)/average_over)*momentum+1.0/average_over*updated_momentum;
}

void Cell::update_energy(double updated_energy) {
    energy = ((average_over-1)/average_over)*energy+1.0/average_over*updated_energy;
}

void Cell::update_temperature(double updated_temperature) {
    temperature = ((average_over-1)/average_over)*temperature+1.0/average_over*updated_temperature;
}

void Cell::update_volume() {
    // Update the effective cell volume. A cell may contain 50% of solid material
    volume = system->width*system->height/(system->settings->cells_x*system->settings->cells_y)*(float)pixels/total_pixels;
}

int Cell::prepare() {
    //* Determine number of candidate collision pairs to be selected in this cell
    double select = system->coeff*particles*(particles-1)*vr_max;

    collision_pairs = round(select);      // Number of pairs to be selected
    return collision_pairs;
}

int Cell::collide(Random *rnd) {
	//* Skip cells with only one particle
    if( particles < 1 ) return 0;  // Skip to the next cell

    vector<Molecule*>&molecules = system->molecules;

    double crm = vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
    int isel, collisions = 0, k, kk, ip1, ip2;
    double cr;

    Molecule *molecule1, *molecule2;

    for( isel=0; isel<collision_pairs; isel++ ) {
		//* Pick two particles at random out of this cell
        k = (int)(rnd->nextDouble()*particles);
        kk = ((int)(k+1+rnd->nextDouble()*(particles-1))) % particles;
        ip1 = particle_indices[k];
        ip2 = particle_indices[kk];

        molecule1 = molecules[ip1];
        molecule2 = molecules[ip2];

		//* Calculate pair's relative speed

        // cr = norm(molecule1->v-molecule2->v,2);
        cr = sqrt(molecule1->v(0)*molecule2->v(0)+molecule1->v(1)*molecule2->v(1)+molecule1->v(2)*molecule2->v(2));

        if( cr > crm ) {         // If relative speed larger than crm,
            crm = cr;            // then reset crm to larger value
        }

		//* Accept or reject candidate pair according to relative speed
        if( cr > rnd->nextDouble()*vr_max ) {
			//* If pair accepted, select post-collision velocities
			collisions++;
            molecule1->collide_with(molecule2,rnd, cr);
		} // Loop over pairs
	}
	
    vr_max = crm;

	return collisions;
}
