#include "cell.h"
#include <math.h>
#include <time.h>
#include <molecule.h>

using namespace std;

Cell::Cell(System *_system) {
    system = _system;
    vr_max = 0;
    particles = 0;
    momentum_time_steps = 0;
    T = system->T;
    momentum = zeros<vec>(3,1);
    momentum_change = zeros<vec>(3,1);
    energy = 0;
    density = 0;
    particle_capacity = 100;
    particle_indices = new unsigned int[particle_capacity];

    volume = system->width*system->height/(system->cells_x*system->cells_y);
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

void Cell::sampleStatistics() {
	// Calculate energy
    energy = 0;
    density = 0;
    vec this_momentum = zeros<vec>(3,1);

    pressure = 0;
    f_sum = 0;
    T = 0;

	Molecule *molecule;
    int particleIndex;
    int atoms_per_molecule;
    for(int i=0;i<particles;i++) {
        particleIndex = particle_indices[i];
        molecule = system->molecules[particleIndex];
        if(!molecule->active) continue;

        energy   += 0.5*dot(molecule->v,molecule->v);
        T += 2*energy/(3*particles);
        if(particles > 10)
            this_momentum += molecule->atoms*molecule->v/particles;
        density  += molecule->atoms;

        atoms_per_molecule = molecule->atoms;
	}

    momentum = (19.0/20)*momentum+1.0/20*this_momentum;

    if(this->particles) {
        f_sum = -2*T;
        // pressure = 1/volume*(particles*T + 0.5*f_sum);
        pressure = particles/volume*T;
    }
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

    Molecule **molecules = system->molecules;

    double crm = vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
    int isel, collisions = 0, k, kk, ip1, ip2;
    double cos_th, sin_th,cr;

	Molecule *molecule1, *molecule2;
	vec vcm = zeros<vec>(3,1);
	vec vrel = zeros<vec>(3,1);

    for( isel=0; isel<collision_pairs; isel++ ) {
		//* Pick two particles at random out of this cell
        k = (int)(rnd->nextDouble()*particles);
        kk = ((int)(k+1+rnd->nextDouble()*(particles-1))) % particles;
        ip1 = particle_indices[k];
        ip2 = particle_indices[kk];

		molecule1 = molecules[ip1];
		molecule2 = molecules[ip2];

		//* Calculate pair's relative speed
        cr = norm(molecule1->v-molecule2->v,2);

		if( cr > crm )         // If relative speed larger than crm,
		crm = cr;            // then reset crm to larger value

		//* Accept or reject candidate pair according to relative speed
        if( cr > rnd->nextDouble()*vr_max ) {
			//* If pair accepted, select post-collision velocities
			collisions++;
            momentum_change += molecule1->collide_with(molecule2,rnd, cr);
		} // Loop over pairs
	}
	
    vr_max = crm;

	return collisions;
}
