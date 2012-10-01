#include "Wall.h"

Wall::Wall(double y, double T, bool upper, double v_x) {
	this->y = y;
	this->T = T;
	this->v_x = v_x;
	this->upper = upper;
}

bool Wall::isMoleculeOutside(Molecule *molecule) {
	return this->upper ? (molecule->r(1) > this->y) : (molecule->r(1) < this->y);
}

double Wall::timeUntilCollision(double y_old, double v_y) {
	return (this->y-y_old)/v_y;
}