#include <cstddef>
#include <cell.h>
#include <math.h>
#include <time.h>
#include <molecule.h>
#include <system.h>

using namespace std;

Cell::Cell(System *_system) {
    system = _system;
    vr_max = 0;
    pixels = 0;
    molecules.reserve(100);
    num_molecules = 0;
    total_pixels = 0;
    dummy_cell = NULL;
    test_value = 1338;
}

bool Cell::cmp(Cell *c1, Cell *c2) {
    return c1->collision_pairs < c2->collision_pairs;
}

void Cell::update_volume() {
    // Update the effective cell volume. A cell may contain 50% of solid material
    volume = system->volume/(system->cells_x*system->cells_y*system->cells_z)*(float)pixels/total_pixels;
    collision_coefficient = 0.5*system->eff_num*M_PI*system->diam*system->diam*system->dt/volume;
}

int Cell::prepare() {
    //* Determine number of candidate collision pairs to be selected in this cell
    double select = collision_coefficient*num_molecules*(num_molecules-1)*vr_max;

    collision_pairs = round(select);      // Number of pairs to be selected

    return collision_pairs;
}

int Cell::collide(Random *rnd) {

    //* Skip cells with only one particle

    if( num_molecules < 2 ) return 0;  // Skip to the next cell

    double crm = vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
    int isel, collisions = 0, ip1, ip2;
    double cr;

    Molecule *molecule1, *molecule2;

    for( isel=0; isel<collision_pairs; isel++ ) {
		//* Pick two particles at random out of this cell
        ip1 = (int)(rnd->nextDouble()*num_molecules);
        ip2 = ((int)(ip1+1+rnd->nextDouble()*(num_molecules-1))) % num_molecules;

        molecule1 = molecules[ip1];
        molecule2 = molecules[ip2];

		//* Calculate pair's relative speed
        cr = sqrt(pow(molecule1->v[0] - molecule2->v[0],2) + pow(molecule1->v[1] - molecule2->v[1],2) + pow(molecule1->v[2] - molecule2->v[2],2));

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

void Cell::add_molecule(Molecule *m) {
    molecules.push_back(m);
    m->index_in_cell = num_molecules;
    m->cell_index = index;
    num_molecules++;
}

void Cell::remove_molecule(Molecule *m) {
    if(num_molecules>1) {
        // Move the last molecule over here
        molecules[m->index_in_cell] = molecules[num_molecules-1];
        molecules[m->index_in_cell]->index_in_cell = m->index_in_cell;
        molecules.erase(molecules.begin()+num_molecules-1);

    } else {
        molecules.erase(molecules.begin());
    }

    num_molecules--;
}
