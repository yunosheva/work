#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
#include <omp.h>
#include "Header.h"
using namespace std;

int main() {
	system("mkdir VTK");
	LatticeBoltzmannComponents Solver;
	Solver.initialize();

	for (int t = 0; t < 10001; t++) {
		Solver.TimeStep(t);
	}
	return 0;
}