#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <chrono>
using namespace std;

const double R = 8.31446;
double delta_t = 1e-9;
double h = 1e-6;
int N_x = 30;
int N_y = 30;
int N_z = 1;
int FullTime = 20000;
int RelaxTime = 1000;
double teta = 1. / 3. * h * h / delta_t / delta_t;
double tau = 1.;
double A = -0.58;
double T = 350;
double g = 0;
double Full_rho1 = 0;
double Full_rho2 = 0;
double Full_rho3 = 0;
vector<double> Full_velocity_x = { 0, 0, 0 };
vector<double> Full_velocity_y = { 0, 0, 0 };
vector<double> Full_velocity_z = { 0, 0, 0 };
vector<double> omega = { 0.01142, 0.2514, 0.0979 }; 
double wetWalls = 1;
double wetObst = 1.;
double sumG = 0;
int volumeObs = 0;
int numberComponent = 2;
vector<double> Tcr = { 190.564, 469.65, 305.51};
vector<double> p_cr = { 4.5992 * 1e6, 3.3675 * 1e6, 4.8711 * 1e6 };
vector<double> mu = { 0.016, 0.07215, 0.03007 };
vector<double> rho_cr = { 162.66, 232, 200};
vector<double> s = { -0.154 , -0.04183, 0.1002 };
double a_0 = 0.4572793, b_0 = 0.07780669;
vector<vector<double>> k = { {0, 0.03, 0.005}, { 0.03, 0, 0.01}, {0.005, 0.01, 0} };
vector<double> gamma = { 0.432, 1., 0.5 };
vector<double> rho_min = { 100, 100, 100 };
vector<double> rho_max = { -100, -100, -100 };
double percent1 = 0.2;
double percent2 = 0.6;
double rho_mix = 280;
double rho_mix_max, rho_mix_min;
double percent = 0.3;

vector<double> fraction_vapor(3);
vector<double> fraction_fluid(3);
vector<double> rho_vapor(3);
vector<double> rho_fluid(3);
double vol_vapor, vol_fluid;



/* равновесные функции распределения, sp - скалярное произведение, u2 - вектор скорости в квадрате */
vector<double> mixture(double per1, double per2) {
	double rho1, rho2, rho3;
	double a1 = mu[0] * mu[1] * per1 / (per1 * mu[0] + (1 - per1) * mu[1]);
	double a2 = (1 - per2) / mu[1] + per2 / mu[2] + a1 * ((1 - per2) / mu[1] / mu[2] +
		per2 / mu[0] / mu[2] - (1 - per2) / mu[1] / mu[1] - per2 / mu[0] / mu[1]);
	rho3 = 1 / a2 * rho_mix * ((1 - per2) / mu[1] - a1 * ((1 - per2) / mu[1] / mu[1]
		+ per2 / mu[0] / mu[1]));
	rho1 = a1 * (rho_mix / mu[1] + rho3 / mu[2] - rho3 / mu[1]);
	rho2 = rho_mix - rho1 - rho3;
	return { rho1, rho2, rho3 };
}
/* равновесные функции распределения, sp - скалярное произведение, u2 - вектор скорости в квадрате */
double F_e1(double sp, double u2, double w, double rho) {
	return w * rho* (1 + sp / teta + sp * sp / (2. * teta * teta) - u2 / 2. / teta);
}

inline double scalar_product(vector<double> vec1, double vec2, double vec3, double vec4) {
	return vec1[0] * vec2 + vec1[1] * vec3 + vec1[2] * vec4;
}

inline double squaring(double vec1, double vec2, double vec3) {
	return vec1 * vec1 + vec2 * vec2 + vec3 * vec3;
}

double F_e(vector<double> vec1, double vec2, double vec3, double vec4, double w, double rho) {
	return w * rho* (1 + scalar_product(vec1, vec2, vec3, vec4) / (teta)+scalar_product(vec1, vec2, vec3, vec4) * scalar_product(vec1, vec2, vec3, vec4) / (2. * teta * teta) - squaring(vec2, vec3, vec4) / 2. / (teta));
}

/* решеточное уравнение Больцмана */
double F(double f, double f_eq, double f_eq1, double f_eq2) {
	return f + (f_eq - f) / tau + (f_eq1 - f_eq2);
}

/* уравнение состояния Пенга-Робинсона */
double a(const double temperature, double omega) {
	double m = 0.382144 + 1.476905 * omega - 0.134488 * omega * omega;
	double a = pow((1 + m * (1 - sqrt(temperature))), 2);
	return a;
}

double PressurePengRobinson(double rho, const double temperature, double omega) {
	double pressure = 1 / 0.307 * (temperature / (1. / rho - 0.253) -
		1.487 * a(temperature, omega) / (1. / rho / rho + 2 * 0.253 / rho - 0.253 * 0.253));
	return pressure;
}

double PressurePengRobinsonMultyComponent(double rho1, double rho2/*, double rho3*/) {
	double B = 0, A = 0, D = 0, S = 0;
	vector<double> Rho = { rho1, rho2/*, rho3*/};
	for (int numComp = 0; numComp < numberComponent; numComp++) {
		double b = b_0 * R * Tcr[numComp] / p_cr[numComp];
		D += R * Rho[numComp] / mu[numComp];
		B += b * Rho[numComp] / mu[numComp];
		S += Rho[numComp] / mu[numComp] * s[numComp] * b;

		for (int num2 = 0; num2 < numberComponent; num2++) {
			A += Rho[num2] / mu[num2] * Rho[numComp] / mu[numComp] * a_0 * R * R * Tcr[num2] * Tcr[numComp] *
				sqrt(a(T / Tcr[numComp], omega[numComp]) * a(T / Tcr[num2], omega[num2]) / p_cr[numComp] / p_cr[num2]) * (1 - k[numComp][num2]);
		}
	}

	double pressure = D * T / (1 + S - B) - A / (pow((1 + S), 2) + 2 * B * (1 + S) - B * B);
	return pressure;
}

/* initial velocity and change velocity*/
vector<vector<vector<vector<double>>>> ux(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
vector<vector<vector<vector<double>>>> uy(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
vector<vector<vector<vector<double>>>> uz(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
vector<vector<vector<vector<double>>>> dux(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
vector<vector<vector<vector<double>>>> duy(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
vector<vector<vector<vector<double>>>> duz(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
vector<vector<vector<double>>> ux_all(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));
vector<vector<vector<double>>> uy_all(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));
vector<vector<vector<double>>> uz_all(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));




/* initial density, mask and "effective" density*/
vector<vector<vector<vector<double>>>> rho(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
vector<vector<vector<int>>> mask(N_x + 2, vector<vector<int>>(N_y + 2, vector<int>(N_z + 2)));
vector<vector<vector<double>>> Fi(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));
vector<vector<vector<double>>> pressure (N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));
vector<vector<vector<double>>> gamma_rho(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));
vector<vector<vector<double>>> Sum_rho(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));
vector<vector<vector<double>>> g_gravity(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)));

void SaveVTKFile(int tStep)
{
	stringstream fname;
	fname << "VTK/adv_";
	if (tStep < 10) fname << "0";
	if (tStep < 100) fname << "0";
	if (tStep < 1000) fname << "0";
	if (tStep < 10000) fname << "0";
	if (tStep < 100000) fname << "0";
	if (tStep < 1000000) fname << "0";
	if (tStep < 10000000) fname << "0";
	fname << tStep << ".vtk";
	ofstream vtk_file(fname.str().c_str());
	vtk_file << "# vtk DataFile Version 3.0\n";
	vtk_file << "Immiscible displacement\n";
	vtk_file << "ASCII\n";
	vtk_file << "DATASET RECTILINEAR_GRID\nDIMENSIONS " << N_x << " " << N_y << " " << N_z;
	vtk_file << "X_COORDINATES " << N_x << " double\n";
	for (int i = 0; i < N_x; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "Y_COORDINATES " << N_y << " double\n";
	for (int i = 0; i < N_y; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "Z_COORDINATES " << N_z << " double\n";
	for (int i = 0; i < N_z; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "POINT_DATA " << N_x * N_y * N_z << endl;

	vtk_file << "SCALARS rho1 double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int l = 1; l < N_z + 1; l++)
		for (int j = 1; j < N_y + 1; j++)
			for (int i = 1; i < N_x + 1; i++)
				vtk_file << rho[0][i][j][l] << " ";
	vtk_file << endl;

	vtk_file << "VECTORS uflow1 double\n";
	for (int l = 1; l < N_z + 1; l++)
		for (int j = 1; j < N_y + 1; j++)
			for (int i = 1; i < N_x + 1; i++)
				vtk_file << ux[0][i][j][l] + g / 2 << " " << uy[0][i][j][l] << " " << uz[0][i][j][l] << " ";
	vtk_file << endl;

	vtk_file << "SCALARS rho_sum double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int l = 1; l < N_z + 1; l++)
		for (int j = 1; j < N_y + 1; j++)
			for (int i = 1; i < N_x + 1; i++)
				vtk_file << rho[1][i][j][l] + rho[0][i][j][l] /*+ rho[2][i][j][l]*/ << " ";
	vtk_file << endl;

	if (numberComponent >= 2) {
		vtk_file << "VECTORS uflow2 double\n";
		for (int l = 1; l < N_z + 1; l++)
			for (int j = 1; j < N_y + 1; j++)
				for (int i = 1; i < N_x + 1; i++)
					vtk_file << ux[1][i][j][l] + g / 2 << " " << uy[1][i][j][l] << " " << uz[1][i][j][l] << " ";
		vtk_file << endl;

		vtk_file << "SCALARS rho2 double 1\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int l = 1; l < N_z + 1; l++)
			for (int j = 1; j < N_y + 1; j++)
				for (int i = 1; i < N_x + 1; i++)
					vtk_file << rho[1][i][j][l] << " ";
		vtk_file << endl;
	}

	if (numberComponent == 3){
		vtk_file << "SCALARS rho3 double 1\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int l = 1; l < N_z + 1; l++)
			for (int j = 1; j < N_y + 1; j++)
				for (int i = 1; i < N_x + 1; i++)
					vtk_file << rho[2][i][j][l] << " ";
		vtk_file << endl;

		vtk_file << "VECTORS uflow3 double\n";
		for (int l = 1; l < N_z + 1; l++)
			for (int j = 1; j < N_y + 1; j++)
				for (int i = 1; i < N_x + 1; i++)
					vtk_file << ux[2][i][j][l] + g / 2 << " " << uy[2][i][j][l] << " " << uz[2][i][j][l] << " ";
		vtk_file << endl;}

	vtk_file << "SCALARS mask double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int l = 1; l < N_z + 1; l++)
		for (int j = 1; j < N_y + 1; j++)
			for (int i = 1; i < N_x + 1; i++)
				vtk_file << mask[i][j][l] << " ";
	vtk_file << endl;

	vtk_file << "SCALARS pressure double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int l = 1; l < N_z + 1; l++)
		for (int j = 1; j < N_y + 1; j++)
			for (int i = 1; i < N_x + 1; i++)
				vtk_file << pressure[i][j][l] << " ";
	vtk_file << endl;

	vtk_file.close();

	cout << endl << "File " << fname.str() << " written" << endl << endl;
}

int main() {
	system("mkdir VTK");

	/* initial density */
	/*for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				if (mask[i][j][l] == 0) {
					rho[0][i][j][l] = 16.5308;
					rho[1][i][j][l] = 42.0405;
					rho[2][i][j][l] = 26.3311;
				}
			}
		}
	}
	for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				if (sqrt((i - 25) * (i - 25) + (j - 25) * (j - 25) + (l) * (l)) <= 3){
					rho[0][i][j][l] = 17.109;
					rho[1][i][j][l] = 416.75;
					rho[2][i][j][l] = 34.8244;
				}
			}
		}
	}

	/*for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				if (mask[i][j][l] == 0) {
					rho[0][i][j][l] = 136;
					rho[1][i][j][l] = 79;
					rho[2][i][j][l] = 96;
				}
			}
		}
	}

	for (int j = 1; j < N_y + 1; j++) {
		for (int l = 1; l < N_z + 1; l++) {
			rho[0][0][j][l] = 141;
			rho[1][0][j][l] = 84;
			rho[2][0][j][l] = 101;
			rho[0][N_x + 1][j][l] = 116;
			rho[1][N_x + 1][j][l] = 59;
			rho[2][N_x + 1][j][l] = 76;
		}
	}*/
	/*for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				if (mask[i][j][l] == 0) {
					rho[0][i][j][l] = 120;
					rho[1][i][j][l] = 70;
					rho[2][i][j][l] = 84;
				}
			}
		}
	}

	for (int j = 1; j < N_y + 1; j++) {
		for (int l = 1; l < N_z + 1; l++) {
			rho[0][0][j][l] = 141;
			rho[1][0][j][l] = 84;
			rho[2][0][j][l] = 101;
			rho[0][N_x + 1][j][l] = 116;
			rho[1][N_x + 1][j][l] = 59;
			rho[2][N_x + 1][j][l] = 76;
		}
	}*/

	/*vector<double> temp(3, 0.0);
	temp = mixture(percent1, percent2);
	cout << temp[0] << endl;
	std::cout << temp[1] << std::endl;
	std::cout << temp[2] << std::endl;
	for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				
				rho[0][i][j][l] = temp[0];
				rho[1][i][j][l] = temp[1] + 2 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
				rho[2][i][j][l] = temp[2];
			}
		}
	}*/
	/*for (int j = 1; j < N_y + 1; j++)
		for (int i = 1; i < N_x + 1; i++)
			for (int l = 1; l < N_z + 1; l++)
				rho[0][i][j][l] = percent * rho_mix * mu[0] / ((1 - percent) * mu[1] + percent * mu[0]);

	for (int j = 1; j < N_y + 1; j++)
		for (int i = 1; i < N_x + 1; i++)
			for (int l = 1; l < N_z + 1; l++)
				rho[1][i][j][l] = rho_mix - rho[0][i][j][l] + 5 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
*/
	for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x / 3; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				rho[0][i][j][l] = 16.1404;
				if (numberComponent >= 2) {
					rho[1][i][j][l] = 52.1604;
				}
				if (numberComponent >= 3) {
					rho[2][i][j][l] = 26.2063;
				}
			}
		}
	}

	for (int j = 1; j < N_y + 1; j++) {
		for (int i = N_x / 3; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				rho[0][i][j][l] = 18.4673;
				if (numberComponent >= 2) {
					rho[1][i][j][l] = 398.393;
				}
				if (numberComponent >= 3) {
					rho[2][i][j][l] = 34.5578;
				}
			}
		}
	}


	/* the law of conservation of mass */
#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				Full_rho1 += rho[0][i][j][l];
				ux[0][i][j][l] = uy[0][i][j][l] = uz[0][i][j][l] = dux[0][i][j][l] = duy[0][i][j][l] = duz[0][i][j][l] = 0.0;
				Full_rho2 += rho[1][i][j][l];
				ux[1][i][j][l] = uy[1][i][j][l] = uz[1][i][j][l] = dux[1][i][j][l] = duy[1][i][j][l] = duz[1][i][j][l] = 0.0;
				if (numberComponent == 3){
					Full_rho3 += rho[2][i][j][l];
					ux[2][i][j][l] = uy[2][i][j][l] = uz[2][i][j][l] = dux[2][i][j][l] = duy[2][i][j][l] = duz[2][i][j][l] = 0.0;}
			}
		}
	}

	/* mask */

	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			mask[i][j][0] = 1;
			mask[i][j][N_z + 1] = 1;  /* walls*/
		}
	}

	for (int i = 1; i < N_x + 1; i++) {
		for (int l = 1; l < N_z + 1; l++) {
			mask[i][0][l] = 1;
			mask[i][N_y + 1][l] = 1;  /* walls*/
		}
	}
	for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {


				if (sqrt((i - 17) * (i - 17) + (j - 8) * (j -8) + (l) * (l)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 29) * (i - 29) + (j - 7) * (j - 7) + (l) * (l)) <= 2)
				{
					//mask[i][j][l] = 1.;
					volumeObs += 1;
				}
				if (sqrt((i - 35) * (i - 35) + (j) * (j)+(l) * (l)) <= 4)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}
				/*if (sqrt((i - 10) * (i - 10) + (j - 2) * (j - 2) + (l - 3) * (l- 3)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 15) * (i - 15) + (j - 4) * (j - 4) + (l-5) * (l-5)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 22) * (i - 22) + (j - 3) * (j - 3) + (l - 3) * (l - 3)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 15) * (i - 15) + (j - 3) * (j - 3) + (l - 1) * (l - 1)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}*/
			}
		}
	}

	/*for (int iter = 0; iter < 4; iter++) {
		for (int i = 1 + 40 * iter; i < 11 + 40 * iter; i++) {
			for (int j = 50; j < 101; j++) {
				mask[i][j][1] = 1;
			}
		}
	}

	for (int iter = 0; iter < 3; iter++) {
		for (int i = 20 + 40 * iter; i < 31 + 40 * iter; i++) {
			for (int j = 60; j < 111; j++) {
				mask[i][j][1] = 1;
			}
		}
	}

	for (int i = 1; i < 23; i++) {
		for (int j = 1; j < 51; j++) {
			mask[i][j][1] = 1;
		}
	}

	for (int i = 109; i < 131; i++) {
		for (int j = 1; j < 51; j++) {
			mask[i][j][1] = 1;
		}
	}
	for (int i = 28; i < 63; i++) {
		for (int j = 1; j < 51; j++) {
			mask[i][j][1] = 1;
		}
	}
	for (int i = 68; i < 103; i++) {
		for (int j = 1; j < 51; j++) {
			mask[i][j][1] = 1;
		}
	}
	/*std::ifstream in("simpleMask.txt");

	if (in.is_open())
	{
		for (int j = 1; j < N_y + 1; j++) {
			for (int i = 1; i < N_x + 1; i++) {
				in >> mask[i][j][1];
			}
		}
	}

	in.close();
/*#pragma omp parallel for
	for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {

				if (sqrt((i - 30) * (i - 30) + (j - 5) * (j - 5) + (l - 5) * (l - 5)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 25) * (i - 25) + (j - 7) * (j - 7) + (l - 7) * (l - 7)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 35) * (i - 35) + (j - 2) * (j - 2) + (l - 7) * (l - 7)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 38) * (i - 38) + (j - 2) * (j - 2) + (l - 1) * (l - 1)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				/*if (sqrt((i - 53) * (i - 53) + (j - 7) * (j - 7) + (l - 4) * (l - 4)) <= 7)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 78) * (i - 78) + (j - 7) * (j - 7) + (l - 20) * (l - 20)) <= 8)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 67) * (i - 67) + (j - 16) * (j - 16) + (l - 8) * (l - 8)) <= 4)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 85) * (i - 85) + (j - 22) * (j - 22) + (l - 15) * (l - 15)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 43) * (i - 43) + (j - 22) * (j - 22) + (l - 17) * (l - 17)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 43) * (i - 43) + (j - 22) * (j - 22) + (l - 15) * (l - 15)) <= 10)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 45) * (i - 45) + (j - 12) * (j - 12) + (l - 5) * (l - 5)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 48) * (i - 48) + (j - 27) * (j - 27) + (l - 44) * (l - 44)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 84) * (i - 84) + (j - 3) * (j - 3) + (l - 4) * (l - 4)) <= 4)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 88) * (i - 88) + (j - 10) * (j - 10) + (l - 8) * (l - 8)) <= 7)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 86) * (i - 86) + (j - 20) * (j - 20) + (l - 4) * (l - 4)) <= 6)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 35) * (i - 35) + (j - 20) * (j - 20) + (l - 17) * (l - 17)) <= 4)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 40) * (i - 40) + (j - 20) * (j - 20) + (l - 17) * (l - 17)) <= 2)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 70) * (i - 70) + (j - 5) * (j - 5) + (l - 4) * (l - 4)) <= 7)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 89) * (i - 89) + (j - 9) * (j - 9) + (l - 18) * (l - 18)) <= 1)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 90) * (i - 90) + (j - 24) * (j - 24) + (l - 14) * (l - 14)) <= 4)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 72) * (i - 72) + (j - 1) * (j - 1) + (l - 1) * (l - 1)) <= 7)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 77) * (i - 77) + (j - 10) * (j - 10) + (l - 8) * (l - 8)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 43) * (i - 43) + (j - 10) * (j - 10) + (l - 8) * (l - 8)) <= 6)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 55) * (i - 55) + (j - 1) * (j - 1) + (l - 18) * (l - 18)) <= 4)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 72) * (i - 72) + (j - 10) * (j - 10) + (l - 8) * (l - 8)) <= 3)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 47) * (i - 47) + (j - 16) * (j - 16) + (l - 5) * (l - 5)) <= 5)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 39) * (i - 39) + (j - 12) * (j - 12) + (l - 12) * (l - 12)) <= 9)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 55) * (i - 55) + (j - 10) * (j - 10) + (l - 17) * (l - 17)) <= 5)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 78) * (i - 78) + (j - 27) * (j - 27) + (l - 40) * (l - 40)) <= 12)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 48) * (i - 48) + (j - 35) * (j - 35) + (l - 41) * (l - 41)) <= 9)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}

				if (sqrt((i - 60) * (i - 60) + (j - 42) * (j - 42) + (l - 30) * (l - 30)) <= 10)
				{
					mask[i][j][l] = 1.;
					volumeObs += 1;
				}
			}
		}
	}*/

	/* vector of possible velocities*/
	vector<vector<double>> c(19, vector<double>(3));
	vector <int> dx = { 0, 1, -1, 0,  0, 0,  0, 1, -1, -1,  1, 1, -1, -1,  1, 0,  0,  0,  0 };
	vector <int> dy = { 0, 0,  0, 1, -1, 0,  0, 1,  1, -1, -1, 0,  0,  0,  0, 1, -1, -1,  1 };
	vector <int> dz = { 0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1,  1, -1, -1, 1,  1, -1, -1 };
	vector <int> index = {         0, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
	vector <int> opposite_index = { 0, 2, 1, 4, 3, 6, 5, 9, 10, 7,  8, 13, 14, 11, 12, 17, 18, 15, 16 };

	for (int j = 1; j < 19; j++) {
		c[j][0] = h * dx[j] / delta_t;
		c[j][1] = h * dy[j] / delta_t;
		c[j][2] = h * dz[j] / delta_t;
	}
	

	/*coefficients w_k */
	vector<double> w(19);

	w[0] = 1. / 3.;

	for (int i = 1; i < 7; i++) {
		w[i] = 1. / 18.;
	}

	for (int i = 7; i < 19; i++) {
		w[i] = 1. / 36.;
	}

	/* coefficients G_k*/
	vector<double> G(19);

	for (int i = 0; i < 7; i++) {
		G[i] = 1.;
	}

	for (int i = 7; i < 19; i++) {
		G[i] = 1. / 2.;
	}

	/* задаем начальные одночастичные функции f, равновесные функции распределения f_eq */
	vector<vector<vector<vector<vector<double>>>>> f(numberComponent, vector<vector<vector<vector<double>>>>(19, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))); 
	vector<vector<vector<vector<double>>>> f1(19, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))));
	vector<vector<vector<vector<vector<double>>>>> buf(numberComponent, vector<vector<vector<vector<double>>>>( 19, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))));
	
#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				for (int s = 0; s < 19; s++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						f[numComp][s][i][j][l] = F_e(c[s], ux[numComp][i][j][l], uy[numComp][i][j][l], uz[numComp][i][j][l], w[s], rho[numComp][i][j][l]);
					}
				}
			}
		}
	}
	auto begin = chrono::steady_clock::now();
	for (int t = 0; t < FullTime + 1; t++) {

		
		/* movement */

		buf = f;

		for (int numComp = 0; numComp < numberComponent; numComp++) {

			if (t < RelaxTime) {

				for (size_t l = 1; l < N_z + 1; l++) {
					for (size_t j = 1; j < N_y + 1; j++) {
						buf[numComp][1][0][j][l] = buf[numComp][2][1][j][l];
						buf[numComp][2][N_x + 1][j][l] = buf[numComp][1][N_x][j][l];
					};
				};
				for (size_t j = 1; j < N_y + 1; j++) {
					for (size_t l = 1; l < N_z + 1; l++) {
						buf[numComp][7][0][j][l] = buf[numComp][9][1][j + 1][l];
					}
				}
				for (size_t j = 1; j < N_y + 1; j++) {
					for (size_t l = 0; l < N_z + 1; l++) {
						buf[numComp][8][N_x + 1][j][l] = buf[numComp][10][N_x][j + 1][l];
					}
				}
				for (size_t j = 1; j < N_y + 1; j++) {
					for (size_t l = 0; l < N_z + 1; l++) {
						buf[numComp][9][N_x + 1][j][l] = buf[numComp][7][N_x][j - 1][l];
					}
				}
				for (size_t j = 1; j < N_y + 1; j++) {
					for (size_t l = 0; l < N_z + 1; l++) {
						buf[numComp][10][0][j][l] = buf[numComp][8][1][j - 1][l];
					}
				}

				for (size_t l = 1; l < N_z + 1; l++) {
					for (size_t j = 0; j < N_y + 1; j++) {
						buf[numComp][11][0][j][l] = buf[numComp][13][1][j][l + 1];
					}
				}

				for (size_t l = 1; l < N_z + 1; l++) {
					for (size_t j = 0; j < N_y + 1; j++) {
						buf[numComp][12][N_x + 1][j][l] = buf[numComp][14][N_x][j][l + 1];
					}
				}

				for (size_t l = 1; l < N_z + 1; l++) {
					for (size_t j = 0; j < N_y + 1; j++) {
						buf[numComp][13][N_x + 1][j][l] = buf[numComp][11][N_x][j][l - 1];
					}
				}

				for (size_t l = 1; l < N_z + 1; l++) {
					for (size_t j = 0; j < N_y + 1; j++) {
						buf[numComp][14][0][j][l] = buf[numComp][12][1][j][l - 1];
					}
				}

				
			}
			else {
				g = 1e-2;
				buf[numComp][1][0] = buf[numComp][1][1];
				buf[numComp][2][N_x + 1] = buf[numComp][2][N_x];
				buf[numComp][7][0] = buf[numComp][7][1];
				buf[numComp][8][N_x + 1] = buf[numComp][8][N_x];
				buf[numComp][9][N_x + 1] = buf[numComp][9][N_x];
				buf[numComp][10][0] = buf[numComp][10][1];
				buf[numComp][11][0] = buf[numComp][11][1];
				buf[numComp][12][N_x + 1] = buf[numComp][12][N_x];
				buf[numComp][13][N_x + 1] = buf[numComp][13][N_x];
				buf[numComp][14][0] = buf[numComp][14][1];
			}

			/* movement with walls*/

			for (int i = 1; i < N_x + 1; i++) {
				for (int l = 1; l < N_z + 1; l++) {
					buf[numComp][3][i][0][l] = buf[numComp][4][i][1][l];
					buf[numComp][4][i][N_y + 1][l] = buf[numComp][3][i][N_y][l];
				}
			}

			for (int i = 1; i < N_x + 1; i++) {
				for (int j = 1; j < N_y + 1; j++) {
					buf[numComp][5][i][j][0] = buf[numComp][6][i][j][1];
					buf[numComp][6][i][j][N_z + 1] = buf[numComp][5][i][j][N_z];
				}
			}

			for (int i = 0; i < N_x + 1; i++) {
				for (int l = 1; l < N_z + 1; l++) {
					buf[numComp][7][i][0][l] = buf[numComp][9][i + 1][1][l];
					buf[numComp][10][i][N_y + 1][l] = buf[numComp][8][i + 1][N_y][l];
				}
			}
			for (int i = 1; i < N_x + 2; i++) {
				for (int l = 1; l < N_z + 1; l++) {
					buf[numComp][9][i][N_y + 1][l] = buf[numComp][7][i - 1][N_y][l];
					buf[numComp][8][i][0][l] = buf[numComp][10][i - 1][1][l];
				}
			}

			for (int i = 0; i < N_x + 1; i++) {
				for (int j = 1; j < N_y + 1; j++) {
					buf[numComp][11][i][j][0] = buf[numComp][13][i + 1][j][1];
					buf[numComp][14][i][j][N_z + 1] = buf[numComp][12][i + 1][j][N_z];
				}
			}
			for (int i = 1; i < N_x + 2; i++) {
				for (int j = 1; j < N_y + 1; j++) {
					buf[numComp][13][i][j][N_z + 1] = buf[numComp][11][i - 1][j][N_z];
					buf[numComp][12][i][j][0] = buf[numComp][14][i - 1][j][1];
				}
			}

			for (int i = 1; i < N_x + 1; i++) {
				for (int j = 0; j < N_y + 1; j++) {
					buf[numComp][15][i][j][0] = buf[numComp][17][i][j + 1][1];
					buf[numComp][18][i][j][N_z + 1] = buf[numComp][16][i][j + 1][N_z];
				}
			}
			for (int i = 1; i < N_x + 1; i++) {
				for (int j = 1; j < N_y + 2; j++) {
					buf[numComp][17][i][j][N_z + 1] = buf[numComp][15][i][j - 1][N_z];
					buf[numComp][16][i][j][0] = buf[numComp][18][i][j - 1][1];
				}
			}

			for (int i = 1; i < N_x + 1; i++) {
				for (int l = 0; l < N_z + 1; l++) {
					buf[numComp][15][i][0][l] = buf[numComp][17][i][1][l + 1];
					buf[numComp][16][i][N_y + 1][l] = buf[numComp][18][i][N_y][l + 1];
				}
			}
			for (int i = 0; i < N_x + 1; i++) {
				for (int l = 1; l < N_z + 2; l++) {
					buf[numComp][17][i][N_y + 1][l] = buf[numComp][15][i][N_y][l - 1];
					buf[numComp][18][i][0][l] = buf[numComp][16][i][1][l - 1];
				}
			}
		}


	/*	for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				for (int s = 0; s < 19; s++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						buf[numComp][s][0][j][l] = buf[numComp][s][1][j][l] + F_e(c[s], 0, 0, 0, w[s], rho[numComp][0][j][l]) - F_e(c[s], 0, 0, 0, w[s], rho[numComp][1][j][l]);
						buf[numComp][s][N_x + 1][j][l] = buf[numComp][s][N_x][j][l] + F_e(c[s], 0, 0, 0, w[s], rho[numComp][N_x + 1][j][l]) - F_e(c[s], 0, 0, 0, w[s], rho[numComp][N_x][j][l]);
					}
				}
			}
		}
		*/


		
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					for (int s = 0; s < 19; s++) {
						for (int numComp = 0; numComp < numberComponent; numComp++) {
							if (mask[i][j][l] == 1) {
								if (mask[i + dx[s]][j + dy[s]][l + dz[s]] != 1)
									buf[numComp][index[s]][i][j][l] = buf[numComp][opposite_index[s]][i + dx[s]][j + dy[s]][l + dz[s]];
							}
						}
					}
				}
			}
		}

#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					for (int s = 0; s < 19; s++) {
						for (int numComp = 0; numComp < numberComponent; numComp++) {
							if (mask[i][j][l] == 0) {
								f[numComp][s][i][j][l] = buf[numComp][s][i - dx[s]][j - dy[s]][l - dz[s]];
							}
						}
					}
				}
			}
		}

		/* new density */
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						rho[numComp][i][j][l] = f[numComp][0][i][j][l];
						for (int s = 1; s < 19; s++) {
							rho[numComp][i][j][l] += f[numComp][s][i][j][l];
						}
					}
				}
			}
		}

		/* the law of conservation of mass */
		Full_rho1 = 0.;
		Full_rho2 = 0.;
		Full_rho3 = 0.;
		rho_min[0] = 1000.;
		rho_max[0] = -1000.;
		rho_min[1] = 1000.;
		rho_max[1] = -1000.;
		rho_min[2] = 1000.;
		rho_max[2] = -1000.;
		rho_mix_max = -1000;
		rho_mix_min = 1000;

#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					Full_rho1 += rho[0][i][j][l];
					Full_rho2 += rho[1][i][j][l];
					if (numberComponent == 3) {
						Full_rho3 += rho[2][i][j][l];
					}
					gamma_rho[i][j][l] = gamma[0] * rho[0][i][j][l];
					Sum_rho[i][j][l] = rho[0][i][j][l];
					for (int numComp = 1; numComp < numberComponent; numComp++) {
						gamma_rho[i][j][l] += gamma[numComp] * rho[numComp][i][j][l];
						Sum_rho[i][j][l] += rho[numComp][i][j][l];
					}
				}
			}
		}
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					g_gravity[i][j][l] = 0.;
					if (rho_min[0] > rho[0][i][j][l]) rho_min[0] = rho[0][i][j][l];
					if (rho_max[0] < rho[0][i][j][l]) rho_max[0] = rho[0][i][j][l];
					if (rho_min[1] > rho[1][i][j][l]) rho_min[1] = rho[1][i][j][l];
					if (rho_max[1] < rho[1][i][j][l]) rho_max[1] = rho[1][i][j][l];
					if (numberComponent == 3) {
						if (rho_min[2] > rho[2][i][j][l]) rho_min[2] = rho[2][i][j][l];
						if (rho_max[2] < rho[2][i][j][l]) rho_max[2] = rho[2][i][j][l];
					}
					if (rho_mix_min > Sum_rho[i][j][l]) rho_mix_min = Sum_rho[i][j][l];
					if (rho_mix_max < Sum_rho[i][j][l]) rho_mix_max = Sum_rho[i][j][l];
					
				}
			};
		};
		
		for (int numComp = 0; numComp < numberComponent; numComp++) {
			rho_vapor[numComp] = 0.;
			rho_fluid[numComp] = 0.;
			vol_fluid = 0.;
			vol_vapor = 0.;
			for (int i = 1; i < N_x + 1; i++) {
				for (int j = 1; j < N_y + 1; j++) {
					for (int l = 1; l < N_z + 1; l++) {
						
						if ((rho_mix_min + rho_mix_max) / 2. > Sum_rho[i][j][l]) {
							rho_vapor[numComp] += rho[numComp][i][j][l];
							vol_vapor += 1.;
						}
						if ((rho_mix_min + rho_mix_max) / 2. <= Sum_rho[i][j][l]) {
							rho_fluid[numComp] += rho[numComp][i][j][l];
							vol_fluid += 1.;
						}
					}
				}
			}
		}
		
		
		for (int numComp = 0; numComp < numberComponent; numComp++) {
			/*if ((rho_vapor[0] / mu[0] + rho_vapor[1] / mu[1] + rho_vapor[2] / mu[2]) == 0.) {
				fraction_vapor[numComp]
			}*/
			//fraction_vapor[numComp] = rho_vapor[numComp] / mu[numComp] / (rho_vapor[0] / mu[0] + rho_vapor[1] / mu[1] + rho_vapor[2] / mu[2]);
			//fraction_fluid[numComp] = rho_fluid[numComp] / mu[numComp] / (rho_fluid[0] / mu[0] + rho_fluid[1] / mu[1] + rho_fluid[2] / mu[2]);
		}
		
		/* new velocity*/
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						ux[numComp][i][j][l] = f[numComp][1][i][j][l] * c[1][0] / rho[numComp][i][j][l];
						uy[numComp][i][j][l] = f[numComp][1][i][j][l] * c[1][1] / rho[numComp][i][j][l];
						uz[numComp][i][j][l] = f[numComp][1][i][j][l] * c[1][2] / rho[numComp][i][j][l];
						for (int s = 2; s < 19; s++) {
							ux[numComp][i][j][l] += f[numComp][s][i][j][l] * c[s][0] / rho[numComp][i][j][l] ;
							uy[numComp][i][j][l] += f[numComp][s][i][j][l] * c[s][1] / rho[numComp][i][j][l];
							uz[numComp][i][j][l] += f[numComp][s][i][j][l] * c[s][2] / rho[numComp][i][j][l];
						}
					}
				}
			}
		}
		
		/* pressure */
	
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
						pressure[i][j][l] = PressurePengRobinsonMultyComponent(rho[0][i][j][l], rho[1][i][j][l]/*, rho[2][i][j][l]*/);
				}
			}
		}
		
		
		/* "effective" density */
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					Fi[i][j][l] = sqrt(Sum_rho[i][j][l] * teta - pressure[i][j][l]);
					//Fi[i][j][l] = 0.;
				}
			}
		}

		/* "effective" density with walls */

		
		Fi[0][0][0] = wetWalls * Fi[1][1][1];
		Fi[0][N_y + 1][0] = wetWalls * Fi[1][N_y][1];
		Fi[0][0][N_z + 1] = wetWalls * Fi[1][1][N_z];
		Fi[0][N_y + 1][N_z + 1] = wetWalls * Fi[1][N_y][N_z];
		Fi[N_x + 1][0][0] = wetWalls * Fi[N_x][1][1];
		Fi[N_x + 1][N_y + 1][0] = wetWalls * Fi[N_x][N_y][1];
		Fi[N_x + 1][0][N_z + 1] = wetWalls * Fi[N_x][1][N_z];
		Fi[N_x + 1][N_y + 1][N_z + 1] = wetWalls * Fi[N_x][N_y][N_z];

		for (int l = 1; l < N_z + 1; l++) {
			Fi[0][0][l] = wetWalls * Fi[1][1][l];
			Fi[0][N_y + 1][l] = wetWalls * Fi[1][N_y][l];
			Fi[N_x + 1][0][l] = wetWalls * Fi[N_x][1][l];
			Fi[N_x + 1][N_y + 1][l] = wetWalls * Fi[N_x][N_y][l];
		}

		for (int i = 1; i < N_x + 1; i++) {
			Fi[i][0][0] = wetWalls * Fi[i][1][1];
			Fi[i][N_y + 1][0] = wetWalls * Fi[i][N_y][1];
			Fi[i][0][N_z + 1] = wetWalls * Fi[i][1][N_z];
			Fi[i][N_y + 1][N_z + 1] = wetWalls * Fi[i][N_y][N_z];
		}

		for (int j = 1; j < N_y + 1; j++) {
			Fi[0][j][0] = wetWalls * Fi[1][j][1];
			Fi[N_x + 1][j][0] = wetWalls * Fi[N_x][j][1];
			Fi[0][j][N_z + 1] = wetWalls * Fi[1][j][N_z];
			Fi[N_x + 1][j][N_z + 1] = wetWalls * Fi[N_x][j][N_z];
		}

		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				Fi[i][j][0] = wetWalls * Fi[i][j][1];
				Fi[i][j][N_z + 1] = wetWalls * Fi[i][j][N_z];
			}
		}

		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				Fi[i][0][l] = wetWalls * Fi[i][1][l];
				Fi[i][N_y + 1][l] = wetWalls * Fi[i][N_y][l];
			}
		}

		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				if (t < RelaxTime) {
					Fi[0][j][l] = wetWalls * Fi[1][j][l];
					Fi[N_x + 1][j][l] = wetWalls * Fi[N_x][j][l];
				}
				else {
					Fi[0][j][l] = Fi[1][j][l];
					Fi[N_x + 1][j][l] = Fi[N_x][j][l];
				}
			}
		}


#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					if (mask[i][j][l] == 1) {
						sumG = 0;
						Fi[i][j][l] = 0;
						for (int s = 1; s < 19; s++) {
							if (mask[i + dx[s]][j + dy[s]][l + dz[s]] != 1) {
								Fi[i][j][l] += Fi[i + dx[s]][j + dy[s]][l + dz[s]] * G[s];
								sumG += G[s];
							}
						}
						if (sumG != 0)
							Fi[i][j][l] = Fi[i][j][l] * wetObst / sumG;

					}
				}
			}
		}

#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						dux[numComp][i][j][l] = 0.;
						duy[numComp][i][j][l] = 0.;
						duz[numComp][i][j][l] = 0.;
						if (mask[i][j][l] == 0.0) {
							for (size_t s = 1; s < 19; s++) {
								dux[numComp][i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]]  * dx[s] +
									A * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dx[s]) *  delta_t * gamma[numComp] / gamma_rho[i][j][l] / h;
								duy[numComp][i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]]  * dy[s] +
									A * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * Fi[i + dx[s]][j + dy[s]][l + dz[s]]  * dy[s]) * delta_t * gamma[numComp] /gamma_rho[i][j][l]/ h;
								duz[numComp][i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dz[s] +
									A * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * Fi[i + dx[s]][j + dy[s]][l + dz[s]]  * dz[s]) * delta_t * gamma[numComp] / gamma_rho[i][j][l] / h;
							}
						}
					}
				}
			}
		}

		/* the law of conservation of impulse */
		Full_velocity_x[0] = 0;
		Full_velocity_y[0] = 0;
		Full_velocity_z[0] = 0;
		Full_velocity_x[1] = 0;
		Full_velocity_y[1] = 0;
		Full_velocity_z[1] = 0;
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						if (mask[i][j][l] == 0) {
							Full_velocity_x[numComp] += rho[numComp][i][j][l] * (ux[numComp][i][j][l] + (dux[numComp][i][j][l] + g) / 2);
							Full_velocity_y[numComp] += rho[numComp][i][j][l] * (uy[numComp][i][j][l] + (duy[numComp][i][j][l]) * 0.5);
							Full_velocity_z[numComp] += rho[numComp][i][j][l] * (uz[numComp][i][j][l] + duz[numComp][i][j][l] * 0.5);
						}
					}
				}
			}
		}
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					ux_all[i][j][l] = rho[0][i][j][l] * ux[0][i][j][l] / Sum_rho[i][j][l];
					uy_all[i][j][l] = rho[0][i][j][l] * uy[0][i][j][l] / Sum_rho[i][j][l];
					uz_all[i][j][l] = rho[0][i][j][l] * uz[0][i][j][l] / Sum_rho[i][j][l];
					for (int numComp = 1; numComp < numberComponent; numComp++) {
						ux_all[i][j][l] += rho[numComp][i][j][l] * ux[numComp][i][j][l] / Sum_rho[i][j][l];
						uy_all[i][j][l] += rho[numComp][i][j][l] * uy[numComp][i][j][l] / Sum_rho[i][j][l];
						uz_all[i][j][l] += rho[numComp][i][j][l] * uz[numComp][i][j][l] / Sum_rho[i][j][l];

					}
				}
			}
		}

		/* one iteration*/
#pragma omp parallel for
		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				for (int l = 1; l < N_z + 1; l++) {
					for (int s = 0; s < 19; s++) {
						for (int numComp = 0; numComp < numberComponent; numComp++) {
							f[numComp][s][i][j][l] = F(f[numComp][s][i][j][l],
								F_e(c[s], ux_all[i][j][l], uy_all[i][j][l], uz_all[i][j][l], w[s], rho[numComp][i][j][l]),
								F_e(c[s], ux[numComp][i][j][l] + dux[numComp][i][j][l] + g, uy[numComp][i][j][l] + duy[numComp][i][j][l], uz[numComp][i][j][l] + duz[numComp][i][j][l], w[s], rho[numComp][i][j][l]),
								F_e(c[s], ux[numComp][i][j][l] , uy[numComp][i][j][l], uz[numComp][i][j][l], w[s], rho[numComp][i][j][l]));
						}
					}
				};
			};
		};

		if (t % 100 == 0)
		{
			SaveVTKFile(t);
			cout << "Summa mass 1 = " << Full_rho1 << endl;
			cout << "Summa mass 2 = " << Full_rho2 << endl;
			cout << "Summa mass 3 = " << Full_rho3 << endl;
			cout << " Impulse for " << t << " step (x) = " << Full_velocity_x[0] << " and  "  <<Full_velocity_x[1] << endl;
			cout << " Impulse for " << t << " step (y) = " << Full_velocity_y[0] << " and  " << Full_velocity_y[1] << endl;
			cout << " Impulse for " << t << " step (z) = " << Full_velocity_z[0] << " and  " << Full_velocity_z[1] << endl;
			cout << " rho mix min = " << rho_mix_min << " and rho mix max = " << rho_mix_max << endl;
			cout << " rho1 min = " << rho_min[0] << " and rho1 max = " << rho_max[0] << endl;
			cout << " rho2 min = " << rho_min[1] << " and rho2 max = " << rho_max[1] << endl;
			if (numberComponent == 3) {
				cout << " rho3 min = " << rho_min[2] << " and rho3 max = " << rho_max[2] << endl;
			}
			cout << " vapor mixture1 = " << fraction_vapor[0] << " and fluid mixture1 = " << fraction_fluid[0] << endl ;
			cout << " vapor mixture2 = " << fraction_vapor[1] << " and fluid mixture2 = " << fraction_fluid[1] << endl;
			if (numberComponent == 3) {
				cout << " vapor mixture3 = " << fraction_vapor[2] << " and fluid mixture3 = " << fraction_fluid[2] << endl;
			}
			/*cout << g << endl;
			if (rho_max[0] - rho_min[0] < 0.1) {
				auto end = std::chrono::steady_clock::now();

				auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
				std::cout << "The time: " << elapsed_ms.count() << " ms\n";
				exit(1);
			}
			/*cout << " fraction 1 = " << (rho_vapor[0] + rho_fluid[0]) / mu[0] / ((rho_vapor[0] + rho_fluid[0]) / mu[0] + (rho_vapor[1] + rho_fluid[1]) / mu[1] + (rho_vapor[2] + rho_fluid[2]) / mu[2])<< endl ;
			cout << " fraction 2 = " << (rho_vapor[1] + rho_fluid[1]) / mu[1] / ((rho_vapor[0] + rho_fluid[0]) / mu[0] + (rho_vapor[1] + rho_fluid[1]) / mu[1] + (rho_vapor[2] + rho_fluid[2]) / mu[2])<< endl ;
			cout << " fraction 1 = " << (rho_vapor[2] + rho_fluid[2]) / mu[2] / ((rho_vapor[0] + rho_fluid[0]) / mu[0] + (rho_vapor[1] + rho_fluid[1]) / mu[1] + (rho_vapor[2] + rho_fluid[2]) / mu[2])<< endl ;
			cout << "Volume vapor = " << vol_vapor << " and volume fluid = " << vol_fluid << endl;
			cout << " smth = " << (rho_vapor[0] * vol_vapor / mu[0] + rho_vapor[1] * vol_vapor / mu[1] + rho_vapor[2] * vol_vapor / mu[2])/ (rho_fluid[0]* vol_fluid / mu[0] + (vol_fluid * rho_fluid[1]) / mu[1] + (vol_fluid * rho_fluid[2]) / mu[2] + rho_vapor[0] * vol_vapor / mu[0] + rho_vapor[1] * vol_vapor / mu[1] + rho_vapor[2] * vol_vapor / mu[2]) << endl;
			cout << " 1 - smth = " << 1 - (rho_vapor[0] * vol_vapor / mu[0] + rho_vapor[1] * vol_vapor / mu[1] + rho_vapor[2] * vol_vapor / mu[2]) / (rho_fluid[0] * vol_fluid / mu[0] + (vol_fluid * rho_fluid[1]) / mu[1] + (vol_fluid * rho_fluid[2]) / mu[2] + rho_vapor[0] * vol_vapor / mu[0] + rho_vapor[1] * vol_vapor / mu[1] + rho_vapor[2] * vol_vapor / mu[2]) << endl;*/
		}

	}
	return 0;
}