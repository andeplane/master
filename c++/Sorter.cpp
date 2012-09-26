#include "Sorter.h"

inline int calcCellIndex(int i, int j, int k, int cellsPerDimension) {
	return k*cellsPerDimension^2 + j*cellsPerDimension + i;
}

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
	int cellsPerDimension = pow(ncell,1.0/3);
	
	int *jx = new int[N];
	int *jy = new int[N];
	int *jz = new int[N];

	for(int n=0;n<N;n++) {
		jx[n] = 0;
		jy[n] = 0;
		jz[n] = 0;
	}

	Molecule *molecule;
	int i, j, k, cell_index;

	//* Count the number of particles in each cell
	for(int n=0;n<ncell;n++)
		this->cell_n[n] = 0;

	for(int n=0; n<N; n++ ) {
		molecule = this->system->molecules[n];

		i = (int)(molecule->r(0)*cellsPerDimension/L); // Cell index for this particle
		j = (int)(molecule->r(1)*cellsPerDimension/L); // Cell index for this particle
		k = (int)(molecule->r(2)*cellsPerDimension/L); // Cell index for this particle
		i = max(i,0);
		j = max(j,0);
		k = max(k,0);
		

		jx[n] = ( i < ncell ) ? i : cellsPerDimension-1;
		jy[n] = ( j < ncell ) ? j : cellsPerDimension-1;
		jz[n] = ( k < ncell ) ? k : cellsPerDimension-1;

		cell_index = calcCellIndex(i,j,k,cellsPerDimension);
		this->cell_n[cell_index]++;
	}
	
	//* Build index list as cumulative sum of the 
	//  number of particles in each cell
	int m=0;
	for(jcell=0; jcell<ncell; jcell++ ) {
		this->index[jcell] = m;
		m += this->cell_n[jcell];
	}

	//* Build cross-reference list
	int *temp = new int [ncell];	  // Temporary array
	for(jcell=0; jcell<ncell; jcell++ )
		temp[jcell] = 0;
	
	for(int n=0; n<N; n++ )	{
		i = jx[n];
		j = jy[n];
		k = jz[n];
		
		cell_index = calcCellIndex(i,j,k,cellsPerDimension);
		
		int idx = this->index[cell_index] + temp[cell_index];
		
		this->Xref[idx] = n;
		temp[cell_index] += 1;
	}
	

	delete [] jx;
	delete [] jy;
	delete [] jz;
	delete [] temp;
}