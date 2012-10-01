#include "Cell.h"
#include <math.h>
#include "lib.h"
#include <time.h>
#include "omp.h"

using namespace std;

Cell::Cell(System *system) {
	this->system = system;
	this->vr_max = 0;
	this->particles = 0;
	this->momentum = zeros<vec>(3,1);
	this->energy = 0;
	this->density = 0;
	double cellLength = system->L / system->cellsPerDimension;
	this->volume = cellLength*cellLength;
}

void Cell::reset() {
	this->particles = 0;
	this->firstParticleIndex = 0;
}

void Cell::sampleStatistics() {
	// Calculate energy
	this->energy = 0;
	this->density = 0;
	this->momentum.zeros();
	Molecule *molecule;
	for(int i=0;i<this->particles;i++) {
		molecule = this->system->molecules[this->firstParticleIndex+i];
		
		this->energy   += 0.5*molecule->atoms*dot(molecule->v,molecule->v);
		this->momentum += molecule->atoms*molecule->v;
		this->density  += molecule->atoms;
	}

	this->energy /= this->volume;
	this->density /= this->volume;
	this->momentum /= this->volume;
}

int Cell::collide() {
	//* Skip cells with only one particle
	if( this->particles < 1 ) return 0;  // Skip to the next cell

	long *idum = system->idums[omp_get_thread_num()];
	
	Sorter *sorter = this->system->sorter;
	Molecule **molecules = this->system->molecules;

	//* Determine number of candidate collision pairs to be selected in this cell
	double select = this->system->coeff*this->particles*(this->particles-1)*this->vr_max;

	int nsel = round(select);      // Number of pairs to be selected
	double crm = this->vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
	int isel, collisions = 0;
	double cos_th, sin_th;

	Molecule *molecule1, *molecule2;
	vec vcm = zeros<vec>(3,1);
	vec vrel = zeros<vec>(3,1);

	for( isel=0; isel<nsel; isel++ ) {
		//* Pick two particles at random out of this cell
		int k = (int)(ran0(idum)*this->particles);
		int kk = ((int)(k+1+ran0(idum)*(this->particles-1))) % this->particles;

		int ip1 = sorter->Xref[ k+this->firstParticleIndex ];      // First particle index
		int ip2 = sorter->Xref[ kk+this->firstParticleIndex ];     // Second particle index
		molecule1 = molecules[ip1];
		molecule2 = molecules[ip2];

		//* Calculate pair's relative speed
		double cr = norm(molecule1->v-molecule2->v,2);

		if( cr > crm )         // If relative speed larger than crm,
		crm = cr;            // then reset crm to larger value

		//* Accept or reject candidate pair according to relative speed
		if( cr > ran0(idum)*this->vr_max ) {
			//* If pair accepted, select post-collision velocities
			collisions++;

			vcm = 0.5*(molecule1->v+molecule2->v);

			cos_th = 1.0 - 2.0*ran0(idum);      // Cosine and sine of
			sin_th = sqrt(1.0 - cos_th*cos_th); // collision angle theta

			vrel(0) = cr*cos_th;             // Compute post-collision
			vrel(1) = cr*sin_th;    // relative velocity
			
			molecule1->v = vcm + 0.5*vrel;
			molecule2->v = vcm - 0.5*vrel;
		} // Loop over pairs
	}
	
	this->vr_max = crm;

	return collisions;
}