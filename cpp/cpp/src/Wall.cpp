#include "Wall.h"

Wall::Wall(double _y, double _T, bool _upper, double _v_x) {
    y = _y;
    T = _T;
    v_x = _v_x;
    upper = _upper;
}

bool Wall::isMoleculeOutside(Molecule *molecule) {
    return upper ? (molecule->r(1) > y) : (molecule->r(1) < y);
}

double Wall::timeUntilCollision(double y_old, double v_y) {
    return (y-y_old)/v_y;
}
