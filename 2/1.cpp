#include <iostream>
#include "omp.h"
#include <ctime>
#include <cstdlib>
#include <string>
using namespace std;

// ** Const **
const unsigned int matrixN = 5;

void printMatrix(double** Mat) {

	for (int i = 0; i<matrixN; i++) {
		for (int j = 0; j<matrixN; j++) {
			cout << Mat[i][j] <<  " ";
		}
		cout << endl;
	}
}

void initMatrixZeros(double** Mat) {
	for (int i = 0; i<matrixN; i++) {
		for (int j = 0; j<matrixN; j++) {
			Mat[i][j] = 0;
		}
	}
}

void initMatrix(double** Mat) {
	for (int i = 0; i<matrixN; i++) {
		for (int j = 0; j<matrixN; j++) {
			Mat[i][j] = rand() % 10;
		}
	}
}

void matrixMult(double** A, double** B, double** C) {

	for (int i = 0; i<matrixN; i++) {
		for (int j = 0; j<matrixN; j++) {
			for (int k = 0; k<matrixN; k++) {
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

int main() {

	srand(time(NULL));

	double** A;
	double** B;
	double** C;

	A = new double *[matrixN];
	B = new double *[matrixN];
	C = new double *[matrixN];

	for(int i = 0; i < matrixN; i++) {
		A[i] = new double[matrixN];
		B[i] = new double[matrixN];
		C[i] = new double[matrixN];
	}

	initMatrix(A);
	initMatrix(B);
	initMatrixZeros(C);

	matrixMult(A, B, C);

	printMatrix(A);
	cout << endl << " '*' " << endl;
	printMatrix(B);
	cout << endl << " '=' " << endl;
	printMatrix(C);
}
