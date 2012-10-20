#include "CollisionObject.h"
#include "System.h"
#include "Molecule.h"
#include "omp.h"

Box::Box(System *system, vec center, double width, double height, double T)  {
	this->system = system;
	this->center = center;
	this->width = width;
	this->height = height;
	this->T = T;

	this->normals = zeros<vec>(4,1);
	this->normals << 1 << -1 << 1 << -1;
	this->v_indices = zeros<mat>(2,4);
	this->v_indices << 1 << 1 << 0 << 0 << endr
					<< 0 << 0 << 1 << 1 << endr;
	this->active = false;
}

bool Box::collide(Molecule *molecule, double dt) { 
    /*
	vec local = molecule->r-this->center;
	if(abs(local(0)) > this->width/2 || abs(local(1)) > this->height/2) {
		if(!this->active) return true;

		// We did collide
		molecule->r -= molecule->v*dt;
		double minTau=1e1000;
		int wall = -1;

		double tau = ((this->center(0)+this->width/2)-molecule->r(0))/molecule->v(0);
		if(tau>0) { minTau=tau; wall = 0; }

		double tau1 = ((this->center(0)-this->width/2)-molecule->r(0))/molecule->v(0);
		if(tau1>0 && tau1 < minTau) { minTau=tau1; wall = 1; }

		double tau2 = ((this->center(1)+this->height/2)-molecule->r(1))/molecule->v(1);
		if(tau2>0 && tau2 < minTau) { minTau=tau2; wall = 2; }

		double tau3 = ((this->center(1)-this->height/2)-molecule->r(1))/molecule->v(1);
		if(tau3>0 && tau3 < minTau) { minTau=tau3; wall = 3; }

		molecule->r += molecule->v*minTau;

		molecule->v(this->v_indices(0,wall)) = sqrt(3.0/2*this->T)*this->system->rand_gauss(idum);
		molecule->v(this->v_indices(1,wall)) = this->normals(wall)*sqrt(-6.0/2*this->T*log(ran0(idum)));

  		molecule->r += molecule->v*(dt-tau);

  		return true;
	}
*/
	return false;
}
