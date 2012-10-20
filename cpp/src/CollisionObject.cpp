#include "CollisionObject.h"
#include "System.h"
#include "Molecule.h"
#include "omp.h"

Box::Box(System *_system, vec _center, double _width, double _height, double _T)  {
    system = _system;
    center = _center;
    width = _width;
    height = _height;
    T = _T;

    normals = zeros<vec>(4,1);
    normals << 1 << -1 << 1 << -1;
    v_indices = zeros<mat>(2,4);
    v_indices << 1 << 1 << 0 << 0 << endr
					<< 0 << 0 << 1 << 1 << endr;
    active = false;
}

bool Box::collide(Molecule *molecule, double dt) { 
    /*
    vec local = molecule->r-center;
    if(abs(local(0)) > width/2 || abs(local(1)) > height/2) {
        if(!active) return true;

		// We did collide
		molecule->r -= molecule->v*dt;
		double minTau=1e1000;
		int wall = -1;

        double tau = ((center(0)+width/2)-molecule->r(0))/molecule->v(0);
		if(tau>0) { minTau=tau; wall = 0; }

        double tau1 = ((center(0)-width/2)-molecule->r(0))/molecule->v(0);
		if(tau1>0 && tau1 < minTau) { minTau=tau1; wall = 1; }

        double tau2 = ((center(1)+height/2)-molecule->r(1))/molecule->v(1);
		if(tau2>0 && tau2 < minTau) { minTau=tau2; wall = 2; }

        double tau3 = ((center(1)-height/2)-molecule->r(1))/molecule->v(1);
		if(tau3>0 && tau3 < minTau) { minTau=tau3; wall = 3; }

		molecule->r += molecule->v*minTau;

        molecule->v(v_indices(0,wall)) = sqrt(3.0/2*T)*system->rand_gauss(idum);
        molecule->v(v_indices(1,wall)) = normals(wall)*sqrt(-6.0/2*T*log(ran0(idum)));

  		molecule->r += molecule->v*(dt-tau);

  		return true;
	}
*/
	return false;
}
