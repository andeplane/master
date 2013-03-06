#include <cell.h>
#include <math.h>
#include <time.h>
#include <molecule.h>
#include <system.h>

using namespace std;

Cell::Cell(System *_system) {
    system = _system;
    vr_max = 0;
    particles = 0;
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

int lolz = 0;

void Cell::update_volume() {
    // Update the effective cell volume. A cell may contain 50% of solid material
    volume = system->Lx*system->Ly*system->Lz/(system->settings->cells_x*system->settings->cells_y*system->settings->cells_z)*(float)pixels/total_pixels;
    collision_coefficient = 0.5*system->eff_num*M_PI*system->diam*system->diam*system->dt/volume;
}

int Cell::prepare() {
    //* Determine number of candidate collision pairs to be selected in this cell
    double select = collision_coefficient*particles*(particles-1)*vr_max;

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

        cr = sqrt(molecule1->v[0]*molecule2->v[0]+molecule1->v[1]*molecule2->v[1]+molecule1->v[2]*molecule2->v[2]);

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
