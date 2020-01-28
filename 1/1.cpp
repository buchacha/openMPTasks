// *** command: g++ -fopenmp -o 1.out -std=c++11 1.cpp
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
const unsigned int matrixN = 20000;
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
                        Mat[i][j] = pow(-1, rand()%2) * (rand() % 1000);
                }
        }
}

pair<unsigned int, unsigned int> maxMinNoParallel(double** Mat) {

	bool findMinRow = false;
	bool findMaxMin = false;

	double localMin;
	double maxMin;

	unsigned int iHat, jHat, jHatLocal;

        for (int i = 0; i<matrixN; i++) {
		findMinRow = false;
                for (int j = 0; j<matrixN; j++) {
			if (!findMinRow || localMin>Mat[i][j]) {
				localMin = Mat[i][j];
				findMinRow = true;
				jHatLocal = j;
			}
                }
        	if (!findMaxMin || maxMin < localMin) {
			maxMin = localMin;
			iHat = i;
			jHat = jHatLocal;
			findMaxMin = true;
		}
	}
	return make_pair(iHat, jHat);
}

pair<unsigned int, unsigned int> maxMinParallel(double** Mat) {

	bool findMinRow = false;
	bool findMaxMin = false;

	double localMin;
	double maxMin;

	unsigned int iHat, jHat, jHatLocal;

	# pragma omp parallel for num_threads(threadsN) shared(maxMin, iHat, jHat, findMaxMin) private(findMinRow, localMin, jHatLocal)
        for (int i = 0; i<matrixN; i++) {
		findMinRow = false;
                for (int j = 0; j<matrixN; j++) {
			if (!findMinRow || localMin>Mat[i][j]) {
				localMin = Mat[i][j];
				findMinRow = true;
				jHatLocal = j;
			}
                }

		# pragma omp critical
		{
	        	if (!findMaxMin || maxMin<localMin) {
				maxMin = localMin;
				iHat = i;
				jHat = jHatLocal;
				findMaxMin = true;
			}
		}
	}
	return make_pair(iHat, jHat);
}

int main() {

        srand(time(NULL));

        double** A;

        A = new double *[matrixN];

        for(int i = 0; i < matrixN; i++) {
                A[i] = new double[matrixN];
        }


        initMatrix(A);
//	printMatrix(A);
//	cout << endl;

	unsigned int idx, jdx;
	auto begin = chrono::steady_clock::now();
	tie(idx, jdx) = maxMinNoParallel(A);
	auto end = chrono::steady_clock::now();
        cout << "Time difference no parallel = " <<
                chrono::duration_cast<chrono::microseconds>(end - begin).count() <<
                " ms" << endl;

	cout << endl;
	cout << "No parallel: maxMin A(" << idx << ", " << jdx << ") = " << A[idx][jdx] << endl << endl;

	begin = chrono::steady_clock::now();
	tie(idx, jdx) = maxMinParallel(A);


	end = chrono::steady_clock::now();
        cout << "Time difference parallel = " <<
                chrono::duration_cast<chrono::microseconds>(end - begin).count() <<
                " ms" << endl;

	cout << endl;
	cout << "Parallel: maxMin A(" << idx << ", " << jdx << ") = " << A[idx][jdx] << endl << endl;

}

