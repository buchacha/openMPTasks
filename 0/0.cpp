// *** Task #0 (Integral, parallel for, reduction, num_threads)
// *** compile command: g++ -fopenmp -o 0.out -std=c++11 0.cpp

#include "omp.h"
#include <iostream>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace std;

const double pi = 3.14159265358979323846;

int const integralN = 10000000;
int const threadsN = 8;

double a = 0;
double b = 1;

double taskFunction(double x) {
        return 4.0/(1+x*x);
}

double integralNoParallel() {

        double h = (b-a)/integralN;
        double s = 0;
        for (int i = 0; i<integralN; i++) {
                s = s + taskFunction(0+i*h)*h;
        }
        return s;
}

double integralParallel() {

        double h = (b-a)/integralN;
        double s = 0;

        # pragma omp parallel for num_threads(threadsN) reduction(+: s)
        for (int i = 0; i<integralN; i++) {
                s = s + taskFunction(0+i*h)*h;
        }
        return s;
}


int main()
{
        double resNoParallel;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        resNoParallel = integralNoParallel();
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        cout << setprecision(10) << resNoParallel << endl;
        cout << "Time difference no parallel = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;

        double resParallel;
        begin = chrono::steady_clock::now();
        resParallel = integralParallel();
        end = chrono::steady_clock::now();
        cout << setprecision(10) << resParallel << endl;
        cout << "Time difference parallel = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
}
