#ifndef COLLISIONOBJECT_H
#define COLLISIONOBJECT_H
class System;
class Molecule;

#include "Molecule.h"

using namespace arma;
using namespace std;

class CollisionObject {
private:
	
public:
	CollisionObject() { }
	virtual bool collide(Molecule *molecule, double dt) { return false; }
};

class Box : public CollisionObject {
private:

public:
	vec center;
	double width;
	double height;
	double T;
	bool active;
	vec normals;
	mat v_indices;
	System *system;

	Box(System *system, vec center, double width, double height, double T);

	virtual bool collide(Molecule *molecule, double dt);
};

#endif