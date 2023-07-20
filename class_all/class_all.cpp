#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
#include <omp.h>
#include "Header.h"
#include "matrix_MRT.h"
using namespace std;

int main() {
	system("mkdir VTK");
	LatticeBoltzmannComponents Solver;
	Solver.initialize();
	MRT solv;
	solv.initialize();
	
	/*for (int t = 0; t < 10001; t++) {
		Solver.TimeStep(t);
		Solver.collisions();
		Solver.data(t);
	}*/

	for (int t = 0; t < 10001; t++) {
		solv.TimeStep(t);
		solv.collisions();
		solv.data(t);
	}

	return 0;
}