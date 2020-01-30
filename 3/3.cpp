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
#include <fstream>
using namespace std;

const unsigned int numberRange = 10;

void printMatrix(double** Mat,  unsigned int matrixN) {
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

void initMatrixZeros(double** Mat,  unsigned int matrixN) {
        for (int i = 0; i<matrixN; i++) {
                for (int j = 0; j<matrixN; j++) {
                        Mat[i][j] = 0;
                }
        }
}

void initMatrix(double** Mat,  unsigned int matrixN) {
        for (int i = 0; i<matrixN; i++) {
                for (int j = 0; j<matrixN; j++) {
                        Mat[i][j] = pow(-1, rand()%2) * (rand() % numberRange);
                }
        }
}

void initVector(double* Vec,  unsigned int matrixN) {
        for (int i = 0; i<matrixN; i++) {
	        Vec[i] = pow(-1, rand()%2) * (rand() % numberRange);
        }
}

void initVectorZeros(double* Vec,  unsigned int matrixN) {
        for (int i = 0; i<matrixN; i++) {
	        Vec[i] = 0;
        }
}

void printVector(double* Vec,  unsigned int matrixN) {
	cout << "[ ";
        for (int i = 0; i<matrixN; i++) {
	        cout << right << setw(5) << setfill(' ') << Vec[i];
        }
	cout << " ]";
}

void rowMult(double** Mat, double* Vec, double* Res,  unsigned int matrixN) {

        for (int i = 0; i<matrixN; i++) {
                for (int j = 0; j<matrixN; j++) {
			Res[i] += Mat[i][j]*Vec[j];
		}
        }
}

void colMult(double** Mat, double* Vec, double* Res, unsigned int matrixN) {
	for (int j = 0; j<matrixN; j++) {
		for (int i = 0; i<matrixN; i++) {
			Res[i] += Mat[i][j]*Vec[j];
		}
        }
}

void ribbonMult(double** Mat, double* Vec, double* Res,  unsigned int threadsN, unsigned int matrixN) {

	unsigned int ribbonN = matrixN < threadsN ? matrixN : threadsN;
	unsigned int ribbonHeight = matrixN / ribbonN;
	unsigned int lastribbonHeight = matrixN % ribbonN;

	for (int k = 0; k<ribbonN; k++) {
		for (int i = 0; i<ribbonHeight; i++) {
	                for (int j = 0; j<matrixN; j++) {
				int iLocal = k*ribbonHeight + i;
				Res[iLocal] += Mat[iLocal][j]*Vec[j];
			}
        	}
        }

	for (int i = 0; i<lastribbonHeight; i++) {
		for (int j = 0; j<matrixN; j++) {
			int iLocal = ribbonN*ribbonHeight + i;
			Res[iLocal] += Mat[iLocal][j]*Vec[j];
		}
        }
}

void rowMultP(double** Mat, double* Vec, double* Res, unsigned int threadsN,  unsigned int matrixN) {

	# pragma omp parallel for num_threads(threadsN)
        for (int i = 0; i<matrixN; i++) {
                for (int j = 0; j<matrixN; j++) {
			Res[i] += Mat[i][j]*Vec[j]; //may race on Vec[j]
		}
        }
}

// TODO: colMatP2 using reduction on Res[]
void colMultP1(double** Mat, double* Vec, double* Res, double** colMat, unsigned int threadsN,  unsigned int matrixN) {

	# pragma omp parallel for num_threads(threadsN)
	for (int j = 0; j<matrixN; j++) {
		for (int i = 0; i<matrixN; i++) {
			colMat[i][j] += Mat[i][j]*Vec[j];
		}
        }

	# pragma omp parallel for num_threads(threadsN)
	for (int i = 0; i<matrixN; i++) {
		for (int j = 0; j<matrixN; j++) {
			Res[i] += colMat[i][j];
		}
	}
}

void ribbonMultP(double** Mat, double* Vec, double* Res, unsigned int threadsN, unsigned int matrixN) {

	unsigned int ribbonN = matrixN < threadsN ? matrixN : threadsN;
	unsigned int ribbonHeight = matrixN / ribbonN;
 	unsigned int lastRibbonHeight = matrixN % ribbonN;

	# pragma omp parallel for num_threads(threadsN)
	for (int k = 0; k<ribbonN; k++) {
		for (int i = 0; i<ribbonHeight; i++) {
	                for (int j = 0; j<matrixN; j++) {
				int iLocal = k*ribbonHeight + i;
				Res[iLocal] += Mat[iLocal][j]*Vec[j];
			}
        	}
        }

	for (int i = 0; i<lastRibbonHeight; i++) {
		for (int j = 0; j<matrixN; j++) {
			int iLocal = ribbonN*ribbonHeight + i;
			Res[iLocal] += Mat[iLocal][j]*Vec[j];
		}
        }
}

int main() {

	unsigned int matrixN;
	unsigned int threadsN;

        srand(time(NULL));


        double** A;

	double** colMat;

	double* b;
	double* res1;
	double* res2;
	double* res3;
	double* res1P;
	double* res2P;
	double* res3P;

	ofstream file2("mytable2.txt");
	if (file2.is_open()) {
		unsigned int tWid = 15;
		file2 << right << setw(tWid) << setfill(' ') << "matrixN";
		file2 << right << setw(tWid) << setfill(' ') << "thredsN";
		file2 << right << setw(tWid) << setfill(' ') << "row,ms";
		file2 << right << setw(tWid) << setfill(' ') << "col,ms";
		file2 << right << setw(tWid) << setfill(' ') << "block,ms";
		file2 << right << setw(tWid) << setfill(' ') << "row P,ms";
		file2 << right << setw(tWid) << setfill(' ') << "col P,ms";
		file2 << right << setw(tWid) << setfill(' ') << "ribbon P,ms";
		file2 << endl;

		for (matrixN=300; matrixN<1011; matrixN += 50) {
			file2 << "------------------------------------------------------------------------------------------------------------------------" << endl;
			for (threadsN=1; threadsN<=8; threadsN++) {

 				A = new double* [matrixN];
			        for(int i = 0; i < matrixN; i++) {
			                A[i] = new double[matrixN];
			        }

				colMat = new double* [matrixN];
			        for(int i = 0; i < matrixN; i++) {
			                colMat[i] = new double[matrixN];
				}

				b = new double[matrixN];

				res1 = new double[matrixN];
				res2 = new double[matrixN];
				res3 = new double[matrixN];
				res1P = new double[matrixN];
				res2P = new double[matrixN];
				res3P = new double[matrixN];

			        initMatrix(A, matrixN);

				initVector(b, matrixN);

				initVectorZeros(res1, matrixN);
				initVectorZeros(res2, matrixN);
				initVectorZeros(res3, matrixN);
				initVectorZeros(res1P, matrixN);
				initVectorZeros(res2P, matrixN);
				initVectorZeros(res3P, matrixN);

				auto begin = chrono::steady_clock::now();
		                rowMultP(A, b, res1P, threadsN, matrixN);
		                auto end = chrono::steady_clock::now();
				auto t1 =
					chrono::duration_cast<chrono::milliseconds>(end - begin).count();

				begin = chrono::steady_clock::now();
		                colMultP1(A, b, res2P, colMat, threadsN, matrixN);
		                end = chrono::steady_clock::now();
		                auto t2 =
		                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

		                begin = chrono::steady_clock::now();
		                ribbonMultP(A, b, res3P,  threadsN, matrixN);
		                end = chrono::steady_clock::now();
		                auto t3 =
		                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

				begin = chrono::steady_clock::now();
		                rowMult(A, b, res1, matrixN);
		                end = chrono::steady_clock::now();
		                auto t4 =
		                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

				begin = chrono::steady_clock::now();
		                colMult(A, b, res2, matrixN);
		                end = chrono::steady_clock::now();
		                auto t5 =
		                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

				begin = chrono::steady_clock::now();
		                ribbonMult(A, b, res3, threadsN, matrixN);
		                end = chrono::steady_clock::now();
		                auto t6 =
		                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

				file2 << right << setw(tWid) << setfill(' ') << matrixN;
				file2 << right << setw(tWid) << setfill(' ') << threadsN;
				file2 << right << setw(tWid) << setfill(' ') << t4;
				file2 << right << setw(tWid) << setfill(' ') << t5;
				file2 << right << setw(tWid) << setfill(' ') << t6;
				file2 << right << setw(tWid) << setfill(' ') << t1;
				file2 << right << setw(tWid) << setfill(' ') << t2;
				file2 << right << setw(tWid) << setfill(' ') << t3;
				file2 << endl;
			}
		}

		file2.close();
	}
}
