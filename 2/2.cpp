// *** command: g++ -fopenmp -o 2.out -std=c++11 2.cpp
#include <iostream>
#include "omp.h"
#include <ctime>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <chrono>
using namespace std;

// ** Const **
const unsigned int matrixN = 500;
const unsigned int threadsN = 8;

void printMatrix(double** Mat) {

	for (int i = 0; i<matrixN; i++) {
		for (int j = 0; j<matrixN; j++) {
			cout << left << setw(5) << setfill(' ') << Mat[i][j];
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

void matrixNoParallelMult(double** A, double** B, double** C) {

	for (int i = 0; i<matrixN; i++) {
		for (int j = 0; j<matrixN; j++) {
			for (int k = 0; k<matrixN; k++) {
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

void matrixParallelMult(double** A, double** B, double** C) {

	# pragma omp parallel for
	for (int i = 0; i<matrixN; i++) {
//		# pragma omp parallel for - doesn't speed up 
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

	auto begin = chrono::steady_clock::now();
	matrixNoParallelMult(A, B, C);
        auto end = chrono::steady_clock::now();
	cout << "Time difference no parallel = " <<
		chrono::duration_cast<chrono::microseconds>(end - begin).count() <<
		" ms" << endl;
//
//	cout << endl <<  "*** No parallel multiplication ***" << endl;
//	printMatrix(A);
//	cout << endl << " * " << endl << endl;
//	printMatrix(B);
//	cout << endl << " = " << endl << endl;
//	printMatrix(C);

	initMatrixZeros(C);

	begin = chrono::steady_clock::now();
	matrixParallelMult(A, B, C);
	end = chrono::steady_clock::now();
	cout << "Time difference parallel = " <<
		chrono::duration_cast<chrono::microseconds>(end - begin).count() <<
		" ms" << endl;

//	cout << endl << "*** Parallel multiplication ***" << endl;
//	printMatrix(A);
//	cout << endl << " * " << endl << endl;
//	printMatrix(B);
//	cout << endl << " = " << endl << endl;
//	printMatrix(C);
}
