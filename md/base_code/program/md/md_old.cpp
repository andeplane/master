#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "atom.h"

using namespace arma;
using namespace std;

int N = 14;
double rho = 0.8;
double T = 1.0;

void initialize();
void initPositions();
void initVelocities();
void rescaleVelocities();
void printPositionsToFile();
double gasdev();

double **r;
double **v;
double **a;
int    *type;
double L;

int main(int argc, char *argv[]) {
	if(argc>1) N = atoi(argv[1]);

	initialize();
	printPositionsToFile();

	return 0;
}

void printPositionsToFile() {
	ofstream file;
	file.open("pos.xyz");
	file << N << endl;
	file << "stupid text" << endl;

	for(int n=0;n<N;n++) {
		file << (type[n] ? "Ar " : "H ") << r[n][0] << " " << r[n][1] << " " << r[n][2] << endl;
	}

	file.close();
}

void initialize() {
	r = new double*[N];
	v = new double*[N];
	a = new double*[N];
	type = new int[N];
	for(int i=0;i<N;i++) {
		r[i] = new double[3];
		v[i] = new double[3];
		a[i] = new double[3];
	}

	initPositions();
	initVelocities();

}

void initPositions() {
	L = pow(N/rho,1.0/3);

	int M=1;
	while(4*M*M*M < N) ++M;
	double a = L/M;

	double xCell[4] = {0.25, 0.75, 0.75, 0.25};
	double yCell[4] = {0.25, 0.75, 0.25, 0.75};
	double zCell[4] = {0.25, 0.25, 0.75, 0.75};

	int n = 0;
	for(int x = 0; x < M; x++) {
		for(int y = 0; y < M; y++) {
			for(int z = 0; z < M; z++) {
				for(int k = 0; k < 4; k++) {
					if(n<N) {
						r[n][0] = (x+xCell[k]) * a;
						r[n][1] = (y+yCell[k]) * a;
						r[n][2] = (z+zCell[k]) * a;
						type[n] = k>0;

						++n;
					}
				}
			}
		}
	}
}

double gasdev() {
	static bool available = false;
	static double gset;
	double fac, rsq, v1, v2;
	if(!available) {
		do {
			v1 = 2.0*rand() / double(RAND_MAX) - 1.0;
			v2 = 2.0*rand() / double(RAND_MAX) - 1.0;
			rsq = v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		available = true;
		return v2*fac;
	}
	else {
		available = false;
		return gset;
	}
}

void initVelocities() {
	for(int n=0;n<N;n++) 
		for(int i=0;i<3;i++) 
			v[n][i] = gasdev();


	double vCM[3] = {0,0,0};
	for(int n=0;n<N;n++)
		for(int i=0;i<3;i++)
			vCM[i] += v[n][i];

	for(int i=0;i<3;i++)
		vCM[i] /= N;
	for(int n=0;n<N;n++) 
		for(int i=0;i<3;i++) 
			v[n][i] -= vCM[i];


	rescaleVelocities();
}

void rescaleVelocities() {
	double vSqdSum = 0;
	for(int n=0;n<N;n++)
		for(int i=0;i<3;i++)
			vSqdSum += v[n][i]*v[n][i];
	double lambda = sqrt(3*(N-1)*T/vSqdSum);
	for(int n=0;n<N;n++) 
		for(int i=0;i<3;i++)
			v[n][i] *= lambda;
}