#include "Cell.h"
#include <math.h>

using namespace std;

Cell::Cell(System *system) {
	this->system = system;
	this->vr_max = 0;
	this->selxtra = 0;
}