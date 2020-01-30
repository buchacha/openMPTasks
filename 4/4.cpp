// *** command: g++ -fopenmp -o 4.out -std=c++11 4.cpp
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

void initNotEqualMatrixZeros(double** Mat,  unsigned int matrixN1,  unsigned int matrixN2) {
        for (int i = 0; i<matrixN1; i++) {
                for (int j = 0; j<matrixN2; j++) {
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

void ribbonMult(double** Mat, double* Vec, double* Res, unsigned int threadsN, unsigned int matrixN) {

        unsigned int ribbonN = matrixN < threadsN ? matrixN : threadsN;
        unsigned int ribbonHeight = matrixN / ribbonN;
        unsigned int lastRibbonHeight = matrixN % ribbonN;

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

unsigned int getBlockN(const unsigned int threadsN) {
	return floor(sqrt(threadsN));
}

void blockMult(double** Mat, double* Vec, double* Res, double** tempMat, unsigned int threadsN, unsigned int matrixN) {

	unsigned int blockN = getBlockN(threadsN);
	unsigned int blockSize = matrixN/blockN;
	unsigned int tailSize = matrixN%blockN;

	for (int k1=0; k1<blockN; k1++) {
		for (int k2=0; k2<blockN; k2++) {
			for (int i=0; i<blockSize; i++) {
				for (int j=0; j<blockSize; j++) {
					int iLocal = k1*blockSize + i;
					int jLocal = k2*blockSize + j;
					tempMat[iLocal][k2] += Mat[iLocal][jLocal]*Vec[jLocal];
				}
			}
		}
		for (int i=0; i<blockSize; i++) {
			for (int j=0; j<tailSize; j++) {
				int iLocal = k1*blockSize + i;
				int jLocal = blockN*blockSize + j;
				tempMat[iLocal][blockN] += Mat[iLocal][jLocal]*Vec[jLocal];
			}
		}
	}

	for (int k=0; k<blockN; k++) {
		for (int i=0; i<tailSize; i++) {
			for (int j=0; j<blockSize; j++) {
				int iLocal = blockN*blockSize + i;
				int jLocal = k*blockSize + j;
				tempMat[iLocal][k] += Mat[iLocal][jLocal]*Vec[jLocal];
			}
		}
	}

	for (int i=0; i<tailSize; i++) {
		for (int j=0; j<tailSize; j++) {
			int iLocal = blockN*blockSize + i;
			int jLocal = blockN*blockSize + j;
			tempMat[iLocal][blockN] += Mat[iLocal][jLocal]*Vec[jLocal];
		}
	}


	for (int ii=0; ii<matrixN; ii++) {
		for (int k=0; k<blockN+1; k++) {
			Res[ii] += tempMat[ii][k];
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

void blockMultP(double** Mat, double* Vec, double* Res, double** tempMat, unsigned int threadsN, unsigned int matrixN) {

	unsigned int blockN = getBlockN(threadsN);
	unsigned int blockSize = matrixN/blockN;
	unsigned int tailSize = matrixN%blockN;

        # pragma omp parallel for num_threads(blockN*blockN)
	for (int k=0; k<blockN*blockN; k++) {
		int k1 = k / blockN;
		int k2 = k % blockN ;
		for (int i=0; i<blockSize; i++) {
			for (int j=0; j<blockSize; j++) {
				int iLocal = k1*blockSize + i;
				int jLocal = k2*blockSize + j;
				tempMat[iLocal][k2] += Mat[iLocal][jLocal]*Vec[jLocal];
			}
		}

		if (k2 == blockN-1) {
			for (int i=0; i<blockSize; i++) {
				for (int j=0; j<tailSize; j++) {
					int iLocal = k1*blockSize + i;
					int jLocal = blockN*blockSize + j;
					tempMat[iLocal][blockN] += Mat[iLocal][jLocal]*Vec[jLocal];
				}
			}
		}
	}

        # pragma omp parallel for num_threads(blockN)
	for (int k=0; k<blockN; k++) {
		for (int i=0; i<tailSize; i++) {
			for (int j=0; j<blockSize; j++) {
				int iLocal = blockN*blockSize + i;
				int jLocal = k*blockSize + j;
				tempMat[iLocal][k] += Mat[iLocal][jLocal]*Vec[jLocal];
			}
		}
	}

        # pragma omp parallel for num_threads(tailSize)
	for (int i=0; i<tailSize; i++) {
		for (int j=0; j<tailSize; j++) {
			int iLocal = blockN*blockSize + i;
			int jLocal = blockN*blockSize + j;
			tempMat[iLocal][blockN] += Mat[iLocal][jLocal]*Vec[jLocal];
		}
	}

	for (int ii=0; ii<matrixN; ii++) {
		for (int k=0; k<blockN+1; k++) {
			Res[ii] += tempMat[ii][k];
		}
	}
}

int main() {

        unsigned int matrixN;
        unsigned int threadsN;

        srand(time(NULL));


        double** A;

        double** tempMat;

        double* b;
        double* res1;
        double* res2;
        double* res1P;
        double* res2P;


        ofstream file2("mytable2.txt");
        if (file2.is_open()) {
                unsigned int tWid = 15;
                file2 << right << setw(tWid) << setfill(' ') << "matrixN";
                file2 << right << setw(tWid) << setfill(' ') << "thredsN";
                file2 << right << setw(tWid) << setfill(' ') << "ribbon,ms";
                file2 << right << setw(tWid) << setfill(' ') << "block,ms";
                file2 << right << setw(tWid) << setfill(' ') << "ribbonP,ms";
                file2 << right << setw(tWid) << setfill(' ') << "blockP,ms";
                file2 << endl;

                for (matrixN=300; matrixN<1011; matrixN += 50) {
                        file2 << "------------------------------------------------------------------------------------------------------------------------" << endl;
                        for (threadsN=1; threadsN<=9; threadsN++) {

                                A = new double* [matrixN];
                                for(int i = 0; i < matrixN; i++) {
                                        A[i] = new double[matrixN];
                                }

				int blockN = getBlockN(matrixN);
                                tempMat = new double* [matrixN];
                                for(int i = 0; i < matrixN; i++) {
                                        tempMat[i] = new double[blockN+1];
                                }

                                b = new double[matrixN];

                                res1 = new double[matrixN];
                                res2 = new double[matrixN];
                                res1P = new double[matrixN];
                                res2P = new double[matrixN];

                                initMatrix(A, matrixN);
                                initNotEqualMatrixZeros(tempMat, matrixN, blockN);

                                initVector(b, matrixN);

                                initVectorZeros(res1, matrixN);
                                initVectorZeros(res2, matrixN);
                                initVectorZeros(res1P, matrixN);
                                initVectorZeros(res2P, matrixN);

                                auto begin = chrono::steady_clock::now();
                                ribbonMult(A, b, res1, threadsN, matrixN);
                                auto end = chrono::steady_clock::now();
                                auto t1 =
                                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

                                begin = chrono::steady_clock::now();
                                blockMult(A, b, res2, tempMat, threadsN, matrixN);
                                end = chrono::steady_clock::now();
                                auto t2 =
                                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

                                begin = chrono::steady_clock::now();
                                ribbonMultP(A, b, res1P, threadsN, matrixN);
                                end = chrono::steady_clock::now();
                                auto t3 =
                                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

                                initNotEqualMatrixZeros(tempMat, matrixN, blockN);
                                begin = chrono::steady_clock::now();
                                blockMultP(A, b, res2P, tempMat, threadsN, matrixN);
                                end = chrono::steady_clock::now();
                                auto t4 =
                                        chrono::duration_cast<chrono::milliseconds>(end - begin).count();

                                file2 << right << setw(tWid) << setfill(' ') << matrixN;
                                file2 << right << setw(tWid) << setfill(' ') << threadsN;
                                file2 << right << setw(tWid) << setfill(' ') << t1;
                                file2 << right << setw(tWid) << setfill(' ') << t2;
                                file2 << right << setw(tWid) << setfill(' ') << t3;
                                file2 << right << setw(tWid) << setfill(' ') << t4;
                                file2 << endl;
                        }
                }

                file2.close();
        }

//	matrixN = 21;
//	threadsN = 8;
//
//	A = new double* [matrixN];
//        for(int i = 0; i < matrixN; i++) {
//                A[i] = new double[matrixN];
//        }
//
//	int blockN = getBlockN(threadsN);
//        tempMat = new double* [matrixN];
//        for(int i = 0; i < matrixN; i++) {
//                tempMat[i] = new double[blockN+1];
//        }
//
//        b = new double[matrixN];
//
//        res1 = new double[matrixN];
//        res2 = new double[matrixN];
//
//	for (int i=0; i<999; i++) {
//
//	        initMatrix(A, matrixN);
//	        initNotEqualMatrixZeros(tempMat, matrixN, blockN+1);
//
//	        initVector(b, matrixN);
//
//	        initVectorZeros(res1, matrixN);
//	        initVectorZeros(res2, matrixN);
//
//		ribbonMult(A, b, res1, threadsN, matrixN);
//
//		blockMultP(A, b, res2, tempMat, threadsN, matrixN);
//
//		for (int i=0; i<matrixN; i++) {
//			if (res1[i]!=res2[i]) {
//				cout << "Error in blockMult" << endl;
//				printVector(res1, matrixN);
//				cout << endl;
//				printVector(res2, matrixN);
//				break;
//
//			}
//		}
//	}
}
