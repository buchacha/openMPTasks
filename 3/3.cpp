// *** command: g++ -fopenmp -o 3.out -std=c++11 3.cpp
// *** question: why openmp only speed up 2-3 times while there are 8 threads and it must be linear speed up
#include <iostream>
#include "omp.h"
#include <ctime>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <chrono>
#include <tuple>
#include "math.h"
using namespace std;

// ** Const **
const unsigned int matrixN = 20;
const unsigned int threadsN = 8;
const unsigned int numberRange = 10;

void printMatrix(double** Mat) {
	cout << "[[";
        for (int i = 0; i<matrixN; i++) {
                for (int j = 0; j<matrixN; j++) {
                        cout << right << setw(5) << setfill(' ') << Mat[i][j];
                }
		if (i != matrixN-1) {
                	cout << " ]" << endl;
			cout << " [";
		}
		else {
			cout << " ]]" << endl;
        	}
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
                        Mat[i][j] = pow(-1, rand()%2) * (rand() % numberRange);
                }
        }
}

void initVector(double* Vec) {
        for (int i = 0; i<matrixN; i++) {
	        Vec[i] = pow(-1, rand()%2) * (rand() % numberRange);
        }
}

void initVectorZeros(double* Vec) {
        for (int i = 0; i<matrixN; i++) {
	        Vec[i] = 0;
        }
}

void printVector(double* Vec) {
	cout << "[ ";
        for (int i = 0; i<matrixN; i++) {
	        cout << right << setw(5) << setfill(' ') << Vec[i];
        }
	cout << " ]";
}

void rowMult(double** Mat, double* Vec, double* Res) {

        for (int i = 0; i<matrixN; i++) {
                for (int j = 0; j<matrixN; j++) {
			Res[i] += Mat[i][j]*Vec[j];
		}
        }
}

void colMult(double** Mat, double* Vec, double* Res) {
	for (int j = 0; j<matrixN; j++) {  // change 'for loops' by their places i<->j
		for (int i = 0; i<matrixN; i++) {
			Res[i] += Mat[i][j]*Vec[j];
		}
        }
}

void blockMult(double** Mat, double* Vec, double* Res) {

// crop matrix into horizontal blocks

	unsigned int blockN = matrixN < threadsN ? matrixN : threadsN;
	unsigned int blockHeight = matrixN / blockN;
	unsigned int lastBlockHeight = blockHeight + matrixN % blockN;

	for (int k = 0; k<blockN-1; k++) {
		for (int i = 0; i<blockHeight; i++) {
	                for (int j = 0; j<matrixN; j++) {
				int iLocal = k*blockHeight + i;
				Res[iLocal] += Mat[iLocal][j]*Vec[j];
			}
        	}
        }

	for (int i = 0; i<lastBlockHeight; i++) {
		for (int j = 0; j<matrixN; j++) {
			int iLocal = (blockN-1)*blockHeight + i;
			Res[iLocal] += Mat[iLocal][j]*Vec[j];
		}
        }
}

int main() {

        srand(time(NULL));

        double** A;

        A = new double* [matrixN];

        for(int i = 0; i < matrixN; i++) {
                A[i] = new double[matrixN];
        }


        initMatrix(A);

	double* b;
	double* res1;
	double* res2;
	double* res3;

	b = new double[matrixN];
	res1 = new double[matrixN];
	res2 = new double[matrixN];
	res3 = new double[matrixN];

	initVector(b);

	initVectorZeros(res1);
	initVectorZeros(res2);
	initVectorZeros(res3);

	rowMult(A, b, res1);
	colMult(A, b, res2);
	blockMult(A, b, res3);

	for (int i=0; i<matrixN; i++) {
		if (res1[i]!=res2[i] || res2[i]!=res3[i] || res3[i]!=res1[i]) {
			cout << "Error somewhere!" << endl;
		}
	}

	printMatrix(A);
	cout << endl;
	printVector(res1);
	cout << endl;
	printVector(res2);
	cout << endl;
	printVector(res3);
	cout << endl;
//        auto begin = chrono::steady_clock::now();
//        auto end = chrono::steady_clock::now();
//        cout << "Time difference no parallel = " <<
//                chrono::duration_cast<chrono::microseconds>(end - begin).count() <<
//                " ms" << endl;
}

