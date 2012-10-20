#include "Sorter.h"

inline int calcCellIndex(int i, int j, int k, int cellsPerDimension) {
	return k*cellsPerDimension*cellsPerDimension + j*cellsPerDimension + i;
}

Sorter::Sorter(System *system) {
    this->system = system;
    this->Xref = new int[system->N];
    this->cellCount = new int[system->numberOfCells];
    for(int n=0;n<system->numberOfCells;n++) {
    	this->cellCount[n] = 0;	
    }
}

void Sorter::sort() {
	//* Find the cell address for each particle
	int N = this->system->N;
	double L = this->system->L;
	Cell **cells = this->system->cells;

	int numberOfCells = this->system->numberOfCells;
	int cellsPerDimension = this->system->cellsPerDimension;
	
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
	for(int n=0;n<numberOfCells;n++)
		cells[n]->reset();

	int k_max = 0;
	int cell_index_max = 0;
	for(int n=0; n<N; n++ ) {
		molecule = this->system->molecules[n];
		i = (int)((molecule->r(0)/L)*cellsPerDimension);
		j = (int)((molecule->r(1)/L)*cellsPerDimension);
		k = (int)((molecule->r(2)/L)*cellsPerDimension);
		i = min(max(i,0),cellsPerDimension-1); // ensure we are within the limits
		j = min(max(j,0),cellsPerDimension-1);
		k = min(max(k,0),cellsPerDimension-1);

		jx[n] = i;
		jy[n] = j;
		jz[n] = k;

		cell_index = calcCellIndex(i,j,k,cellsPerDimension);

		this->cellCount[cell_index]++;

		cells[cell_index]->particles++;
	}

	//* Build index list as cumulative sum of the 
	//  number of particles in each cell
	int m=0;
	for(int n=0; n<numberOfCells; n++ ) {
		cells[n]->firstParticleIndex = m;
		m += cells[n]->particles;
	}	

	//* Build cross-reference list
	int *temp = new int [numberOfCells];	  // Temporary array
	for(int n=0; n<numberOfCells; n++ )
		temp[n] = 0;
	
	for(int n=0; n<N; n++ )	{
		i = jx[n];
		j = jy[n];
		k = jz[n];
		
		cell_index = calcCellIndex(i,j,k,cellsPerDimension);
		
		int idx = cells[cell_index]->firstParticleIndex + temp[cell_index];
		
		this->Xref[idx] = n;
		temp[cell_index] += 1;
	}
	

	delete [] jx;
	delete [] jy;
	delete [] jz;
	delete [] temp;
}