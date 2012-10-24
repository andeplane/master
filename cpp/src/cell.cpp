#include "Cell.h"
#include <math.h>
#include <time.h>
#include <molecule.h>

using namespace std;

Cell::Cell(System *_system) {
    system = _system;
    vr_max = 0;
    particles = 0;
    momentum = zeros<vec>(3,1);
    momentum_change = zeros<vec>(3,1);
    energy = 0;
    density = 0;
    particle_capacity = 100;
    particle_indices = new unsigned int[particle_capacity];

    volume = system->width*system->height/(system->cells_x*system->cells_y);
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
    momentum.zeros();
    pressure = 0;

	Molecule *molecule;
    int particleIndex;
    int atoms_per_molecule;
    for(int i=0;i<particles;i++) {
        particleIndex = particle_indices[i];
        molecule = system->molecules[particleIndex];

        energy   += 0.5*dot(molecule->v,molecule->v);

        momentum += molecule->atoms*molecule->v;
        density  += molecule->atoms;

        atoms_per_molecule = molecule->atoms;
	}

    if(this->particles) {
        pressure += norm(momentum_change,2)/2;
        pressure += 2*this->particles*energy/3;

        pressure /= this->particles*system->dt;
        momentum_change.zeros();
    }

    energy /= volume;
    density /= volume;
    momentum /= volume;
}

int Cell::collide(Random *rnd) {
	//* Skip cells with only one particle
    if( particles < 1 ) return 0;  // Skip to the next cell

    Sorter *sorter = system->sorter;
    Molecule **molecules = system->molecules;

	//* Determine number of candidate collision pairs to be selected in this cell
    double select = system->coeff*particles*(particles-1)*vr_max;

	int nsel = round(select);      // Number of pairs to be selected
    double crm = vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
	int isel, collisions = 0;
	double cos_th, sin_th;

	Molecule *molecule1, *molecule2;
	vec vcm = zeros<vec>(3,1);
	vec vrel = zeros<vec>(3,1);

	for( isel=0; isel<nsel; isel++ ) {
		//* Pick two particles at random out of this cell
        int k = (int)(rnd->nextDouble()*particles);
        int kk = ((int)(k+1+rnd->nextDouble()*(particles-1))) % particles;
        int ip1 = particle_indices[k];
        int ip2 = particle_indices[kk];
        // int ip1 = sorter->Xref[ k+firstParticleIndex ];      // First particle index
        // int ip2 = sorter->Xref[ kk+firstParticleIndex ];     // Second particle index
		molecule1 = molecules[ip1];
		molecule2 = molecules[ip2];

		//* Calculate pair's relative speed
		double cr = norm(molecule1->v-molecule2->v,2);

		if( cr > crm )         // If relative speed larger than crm,
		crm = cr;            // then reset crm to larger value

		//* Accept or reject candidate pair according to relative speed
        if( cr > rnd->nextDouble()*vr_max ) {
			//* If pair accepted, select post-collision velocities
			collisions++;

			vcm = 0.5*(molecule1->v+molecule2->v);
            double len = norm(vcm,2);

            cos_th = 1.0 - 2.0*rnd->nextDouble();      // Cosine and sine of
			sin_th = sqrt(1.0 - cos_th*cos_th); // collision angle theta

			vrel(0) = cr*cos_th;             // Compute post-collision
			vrel(1) = cr*sin_th;    // relative velocity

            momentum_change += vcm + 0.5*vrel-molecule1->v;
			
			molecule1->v = vcm + 0.5*vrel;
			molecule2->v = vcm - 0.5*vrel;
            vec new_vcm = 0.5*(molecule1->v+molecule2->v);

            double len2 = norm(new_vcm,2);
            if(abs(len2-len) > 1e-8) {
                cout << "WHAT the fuck is going on???" << endl;
            }
		} // Loop over pairs
	}
	
    vr_max = crm;

	return collisions;
}
