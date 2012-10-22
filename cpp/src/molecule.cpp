#include "Molecule.h"
#include <iostream>
#include <math.h>
#include "Wall.h"
#include "omp.h"
#include <CVector.h>
using namespace std;

Molecule::Molecule(System *_system) {
    r = zeros<vec> (3,1);
    v = zeros<vec> (3,1);

    atoms = 1;
    type = 0;
    system = _system;
}

inline int sign(double a) {
	return a >= 0 ? 1 : -1;
}

void Molecule::move(double dt, Random *rnd, int depth) {
    double tau = dt;
    // cout << "i=" << system->world_grid->get_grid_point(r)->i << ", j=" << system->world_grid->get_grid_point(r)->j << endl;

    r += v*dt;
    // cout << "moved, i=" << system->world_grid->get_grid_point(r)->i << ", j=" << system->world_grid->get_grid_point(r)->j << endl;

    GridPoint *point = system->world_grid->get_grid_point(r);

    // We have to calculate time until collision
    if(point->is_wall) {
        // cout << "Is wall: " << point->i << ", " << point->j << endl<< endl<< endl<< endl;
        // cout << "Is boundary? " << point->is_wall_boundary << endl;
        if(!point->is_wall_boundary) {
            int count = 0;
            while(true) {
                point = system->world_grid->get_grid_point(r);
                // cout << "I am now at " << point->i << ", " << point->j << "is boundary: " << point->is_wall_boundary << endl;

                if(point->is_wall_boundary) {
                    dt -= tau;

                    while(system->world_grid->get_grid_point(r)->is_wall) {
                        dt += 0.1*tau;
                        r -= 0.1*v*tau;
                    }
                    // cout<< "I did break here" << endl;
                    break;
                }

                if(point->is_wall) {
                    r -= v*tau;
                    tau /= 2;
                    // cout << "is wall, moved back to: i=" << system->world_grid->get_grid_point(r)->i << ", j=" << system->world_grid->get_grid_point(r)->j << endl;
                }
                else {
                    // cout << "is not wall, moving forward" << endl;
                    dt -= tau;
                }

                if(++count > 100) {
                    cout << "Damn, we screwed up at " << point->i << " , " << point->j << endl;
                    exit(1);
                }

                r += v*tau;
            }
        }
        else {
            // This was actually the boundary
            while(system->world_grid->get_grid_point(r)->is_wall) {
                dt += 0.1*tau;
                r -= 0.1*v*tau;
            }
        }

        double v_normal = sqrt(-6.0/2*point->T*log(rnd->nextDouble()));
        double v_tangent = sqrt(3.0/2*point->T)*rnd->nextGauss();

        v(0) = v_normal*point->normal.x + v_tangent*point->tangent.x;
        v(1) = v_normal*point->normal.y + v_tangent*point->tangent.y;

    }
    else dt = 0;

    r += v*dt;

    if(system->world_grid->get_grid_point(r)->is_wall) {
        GridPoint *p = system->world_grid->get_grid_point(r);

        cout << "OMFG, WE CAN'T WORK LIKE THIS!" << endl;
        cout << "Sorry, I suck at " << p->i << ", " << p->j << endl;
        cout << "Normal vector: " << endl << p->normal << endl;
        cout << "Tangent vector: " << endl << p->tangent << endl;
        exit(1);
    }

    if(depth == 0) {
        r(0) = fmod(r(0)+10*system->width,system->width);
        r(1) = fmod(r(1)+10*system->height,system->height);
    }
}

void Molecule::move_old(double dt, Random *rnd) {
    double y_old = r(1);
    r += v*dt;

    int outsideWallIndex = system->walls[0]->isMoleculeOutside(this) + 2*system->walls[1]->isMoleculeOutside(this);
	double tau = 1;

	if(outsideWallIndex) {
        Wall *wall = system->walls[--outsideWallIndex]; // Decrease from [1,2] to [0,1] first
		int direction = wall->upper ? -1 : 1;

		// We are outside the box
		tau = wall->timeUntilCollision(y_old,v(1));
		
        r += v*(tau-dt);

        v(0) = sqrt(3.0/2*wall->T)*rnd->nextGauss() + wall->v_x;
        v(1) = direction*sqrt(-6.0/2*wall->T*log(rnd->nextDouble()));

        r += v*(dt-tau);
	}

    r(0) = fmod(r(0)+10*system->width,system->width);
    r(1) = fmod(r(1)+10*system->height,system->height);
}
