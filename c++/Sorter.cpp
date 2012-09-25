#include "Sorter.h"

Sorter::Sorter(System *system) {
    this->system = system;
    this->ncell = system->ncell;


    this->cell_n = new int[this->ncell];
    this->index = new int[this->ncell];
    this->Xref = new int[system->N];
}

void Sorter::sort() {

	//* Find the cell address for each particle
	int N = this->system->N;
	double L = this->system->L;
	int jcell;
	int ncell = this->ncell;

	int *jx = new int[N];
	for(int jcell=0;jcell<ncell;jcell++) 
		jx[jcell] = 0;

	Molecule *molecule;
	for(int n=0; n<N; n++ ) {
		molecule = this->system->molecules[n];

		int j = (int)(molecule->r(0)*ncell/L); // Cell index for this particle
		jx[n] = ( j < ncell ) ? j : ncell-1; // Particle n is in cell j
	}

	//* Count the number of particles in each cell
	for(int n=0;n<ncell;n++)
		this->cell_n[n] = 0;
	
	for(int n=0; n<N; n++)
		this->cell_n[ jx[n] ]++; // Increase particle count in this cell
	
	//* Build index list as cumulative sum of the 
	//  number of particles in each cell
	int m=0;
	for(jcell=0; jcell<ncell; jcell++ ) {
		this->index[jcell] = m;
		m += this->cell_n[jcell];
	}

	//* Build cross-reference list
	int *temp;
	temp = new int [ncell];	  // Temporary array
	for(jcell=0; jcell<ncell; jcell++ )
		temp[jcell] = 0;

	for(int n=0; n<N; n++ )	{
		jcell = jx[n]; // Cell index of particle n
		int k = this->index[jcell] + temp[jcell];
		this->Xref[k] = n;
		temp[jcell] += 1;
	}

	delete [] jx;
	delete [] temp;
}