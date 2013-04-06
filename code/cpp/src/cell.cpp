#include <cstddef>
#include <cell.h>
#include <math.h>
#include <time.h>
#include <molecule.h>
#include <system.h>
#include <structs.h>

using namespace std;

Cell::Cell(System *_system) {
    system = _system;
    vr_max = 0;
    pixels = 0;
    num_molecules = 0;
    total_pixels = 0;
    dummy_cell = NULL;
    r         = new double[3*100000];
    v        = new double[3*100000];
    r0 = new double[3*100000];
    atom_moved        = new bool[100000];
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

void Cell::collide_molecules(const int &ip0, const int &ip1, const double &v_rel, Random *rnd) {
    double vcmx  = 0.5*(v[3*ip0+0] + v[3*ip1+0]);
    double vcmy  = 0.5*(v[3*ip0+1] + v[3*ip1+1]);
    double vcmz  = 0.5*(v[3*ip0+2] + v[3*ip1+2]);
    double vx0 = v[3*ip0+0];
    double vy0 = v[3*ip0+1];
    double vz0 = v[3*ip0+2];

    double cos_th = 1.0 - 2.0*rnd->nextDouble();      // Cosine and sine of
    double sin_th = sqrt(1.0 - cos_th*cos_th);        // collision angle theta
    double phi = 2*M_PI*rnd->nextDouble();

    double vrelx = v_rel*cos_th;                   // Compute post-collision relative velocity
    double vrely = v_rel*sin_th*cos(phi);
    double vrelz = v_rel*sin_th*sin(phi);

    v[3*ip0+0] = vcmx + 0.5*vrelx;
    v[3*ip0+1] = vcmy + 0.5*vrely;
    v[3*ip0+2] = vcmz + 0.5*vrelz;

    v[3*ip1+0] = vcmx - 0.5*vrelx;
    v[3*ip1+1] = vcmy - 0.5*vrely;
    v[3*ip1+2] = vcmz - 0.5*vrelz;
}

int Cell::collide(Random *rnd) {

    //* Skip cells with only one particle

    if( num_molecules < 2 ) return 0;  // Skip to the next cell

    double crm = vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
    int isel, collisions = 0, ip0, ip1;
    double v_rel;

    for( isel=0; isel<collision_pairs; isel++ ) {
		//* Pick two particles at random out of this cell
        ip0 = (int)(rnd->nextDouble()*num_molecules);
        ip1 = ((int)(ip1+1+rnd->nextDouble()*(num_molecules-1))) % num_molecules;

        //* Calculate pair's relative speed

        v_rel = sqrt(pow(v[3*ip0+0] - v[3*ip1+0],2) + pow(v[3*ip0+1] - v[3*ip1+1],2) + pow(v[3*ip0+2] - v[3*ip1+2],2));

        if( v_rel > crm ) {         // If relative speed larger than crm,
            crm = v_rel;            // then reset crm to larger value
        }

		//* Accept or reject candidate pair according to relative speed
        if( v_rel > rnd->nextDouble()*vr_max ) {
			//* If pair accepted, select post-collision velocities
			collisions++;
            collide_molecules(ip0,ip1, v_rel, rnd);
		} // Loop over pairs
	}
	
    vr_max = crm;

	return collisions;
}

void Cell::update_molecule_cells(int dimension) {
    for(int n=0;n<num_molecules;n++) {
        if(atom_moved[n]) continue; // Skip atoms that are already moved

        int cell_index = system->thread_control.cell_index_from_position(&r[3*n]);
        DummyCell *dummy_cell = system->thread_control.dummy_cells[cell_index];

        if(cell_index != this->index && (dummy_cell->index_vector[dimension] == this->dummy_cell->index_vector[dimension])) {
            struct Molecule molecule;
            for(int a=0;a<3;a++) {
                molecule.r[a]  = r[3*n+a];
                molecule.v[a]  = v[3*n+a];
                molecule.r0[a] = r0[3*n+a];
            }
            system->thread_control.add_molecule_to_cell(molecule,cell_index);
            atom_moved[n] = true;
        }
    }
}

void Cell::add_molecule(struct Molecule &m) {
    add_molecule(m.r,m.v,m.r0);
}

void Cell::add_molecule(double *r_, double *v_, double *r0_) {
    r[3*num_molecules + 0] = r_[0];
    r[3*num_molecules + 1] = r_[1];
    r[3*num_molecules + 2] = r_[2];

    r0[3*num_molecules + 0] = r0_[0];
    r0[3*num_molecules + 1] = r0_[1];
    r0[3*num_molecules + 2] = r0_[2];

    v[3*num_molecules + 0] = v_[0];
    v[3*num_molecules + 1] = v_[1];
    v[3*num_molecules + 2] = v_[2];
    atom_moved[num_molecules] = false;

    num_molecules++;
}

void Cell::add_molecule(double *r_, double *v_) {
    r[3*num_molecules + 0] = r_[0];
    r[3*num_molecules + 1] = r_[1];
    r[3*num_molecules + 2] = r_[2];

    r0[3*num_molecules + 0] = r_[0];
    r0[3*num_molecules + 1] = r_[1];
    r0[3*num_molecules + 2] = r_[2];

    v[3*num_molecules + 0] = v_[0];
    v[3*num_molecules + 1] = v_[1];
    v[3*num_molecules + 2] = v_[2];
    atom_moved[num_molecules] = false;

    num_molecules++;
}

void Cell::update_molecule_arrays() {
    int remaining_molecules = 0;
    for(int n=0;n<num_molecules;n++) {
        if(!atom_moved[n]) {
            for(int a=0;a<3;a++) {
                r[3*remaining_molecules+a]  = r[3*n+a];
                v[3*remaining_molecules+a]  = v[3*n+a];
                r0[3*remaining_molecules+a] = r0[3*n+a];
            }

            atom_moved[remaining_molecules] = false;
            remaining_molecules++;
        }
    }

    num_molecules = remaining_molecules;
}
