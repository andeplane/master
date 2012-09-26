#include "Cell.h"
#include <math.h>
#include "lib.h"

using namespace std;

const double pi = 3.141592654;

Cell::Cell(System *system) {
	this->system = system;
	this->vr_max = 0;
	this->selxtra = 0;
	this->particlesInCell = 0;
}

void Cell::reset() {
	this->particlesInCell = 0;
}

int Cell::collide() {
	Sorter *sorter = this->system->sorter;
	Molecule **molecules = this->system->molecules;

	//* Skip cells with only one particle
	if( this->particlesInCell < 1 ) return 0;  // Skip to the next cell

	//* Determine number of candidate collision pairs
	//  to be selected in this cell

	double select = this->system->coeff*this->particlesInCell*(this->particlesInCell-1)*this->vr_max + this->selxtra;

	int nsel = (int)(select);      // Number of pairs to be selected
	this->selxtra = select-nsel;  // Carry over any left-over fraction
	double crm = this->vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
	int isel, col = 0;
	Molecule *molecule1, *molecule2;

	for( isel=0; isel<nsel; isel++ ) {

	  //* Pick two particles at random out of this cell
	  int k = (int)(ran0(this->system->idum)*this->particlesInCell);
	  int kk = ((int)(k+1+ran0(this->system->idum)*(this->particlesInCell-1))) % this->particlesInCell;
	  int ip1 = sorter->Xref[ k+sorter->index[this->index] ];      // First particle
	  int ip2 = sorter->Xref[ kk+sorter->index[this->index] ];     // Second particle
	  molecule1 = molecules[ip1];
	  molecule2 = molecules[ip2];

	  //* Calculate pair's relative speed
	  double cr = norm(molecule1->v-molecule2->v,2);

	  if( cr > crm )         // If relative speed larger than crm,
	    crm = cr;            // then reset crm to larger value

	  //* Accept or reject candidate pair according to relative speed
	  if( cr/this->vr_max > ran0(this->system->idum) ) {
	    //* If pair accepted, select post-collision velocities
	    col++;                     // Collision counter
	    vec vcm = zeros<vec>(3,1);
	    vec vrel = zeros<vec>(3,1);

	    int k;
	    vcm = 0.5*(molecule1->v+molecule2->v);

	    double cos_th = 1.0 - 2.0*ran0(this->system->idum);      // Cosine and sine of
	    double sin_th = sqrt(1.0 - cos_th*cos_th); // collision angle theta
	    double phi = 2.0*pi*ran0(this->system->idum);            // Collision angle phi
	    vrel(0) = cr*cos_th;             // Compute post-collision
	    vrel(1) = cr*sin_th*cos(phi);    // relative velocity
	    vrel(2) = cr*sin_th*sin(phi);

	    molecule1->v = vcm + 0.5*vrel;
	    molecule2->v = vcm - 0.5*vrel;
	  } // Loop over pairs
	}

	this->vr_max = crm;
	return col;
}
