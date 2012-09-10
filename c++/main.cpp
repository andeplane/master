#include <iostream>
#include "time.h"
using namespace std;

void matrMult1(double **A, double **B, double **C,int N) {
	for(int i=0;i<N;i++) {
		for(int j=0;j<N;j++) {
			for(int k=0;k<N;k++) {
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

void matrMult2(double **A, double **B, double **C,int N) {
	for(int i=0;i<N;i++) {
		for(int j=0;j<N;j++) {
			for(int k=0;k<N;k++) {
				C[j][i] += A[j][k]*B[k][i];
			}
		}
	}
}

int main(int args, char* argv[]) {
	double **A = (double**)malloc(N*sizeof(double*));
	double **B = (double**)malloc(N*sizeof(double*));
	double **C1 = (double**)malloc(N*sizeof(double*));
	double **C2 = (double**)malloc(N*sizeof(double*));
}