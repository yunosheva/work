#include "Header.h"
#include "matrix_MRT.h"

LatticeBoltzmannComponents::LatticeBoltzmannComponents() :
ux(vector<vector<vector<vector<double>>>>(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))),
uy(vector<vector<vector<vector<double>>>>(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))),
uz(vector<vector<vector<vector<double>>>>(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))),
dux(vector<vector<vector<vector<double>>>>(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))),
duy(vector<vector<vector<vector<double>>>>(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))),
duz(vector<vector<vector<vector<double>>>>(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))),
rho(vector<vector<vector<vector<double>>>>(numberComponent, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2))))),
f(vector<vector<vector<vector<vector<double>>>>>(numberComponent, vector<vector<vector<vector<double>>>>(kMax, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))))),
buf(vector<vector<vector<vector<vector<double>>>>>(numberComponent, vector<vector<vector<vector<double>>>>(kMax, vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))))),
ux_all(vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))),
uy_all(vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))),
uz_all(vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))),
Fi(vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))),
gamma_rho(vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))),
Sum_rho(vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))),
pressure(vector<vector<vector<double>>>(N_x + 2, vector<vector<double>>(N_y + 2, vector<double>(N_z + 2)))),
mask(vector<vector<vector<int>>>(N_x + 2, vector<vector<int>>(N_y + 2, vector<int>(N_z + 2)))),
Full_velocity_x(vector<double>(numberComponent, 0.0)),
Full_velocity_y(vector<double>(numberComponent, 0.0)),
Full_velocity_z(vector<double>(numberComponent, 0.0)),
Full_rho(vector<double>(numberComponent, 0.0)),
rho_max(vector<double>(numberComponent, -1000.0)),
rho_min(vector<double>(numberComponent, 1000.0)),
c(vector<vector<double>>(kMax, vector<double>(3))),
wall_condition(vector<bool>(3)),
X_condition(vector<bool>(3)),
w(vector<double>(kMax)),
G(vector<double>(kMax))
{
	switch (numberComponent) {
	case 1:
	{
		Tcr = { 190.56 };
		p_cr = { 4.5992 * 1e6 };
		mu = { 0.016, 0.07215 };
		omega = { 0.01142};
		rho_cr = { 162.66};
		s = { -0.154  };
		k = { {0} };
		gamma = { 0.432 };
	}
	case 2: // metan, pentan
	{
		Tcr = { 190.564, 469.65 }; 
		p_cr = { 4.5992 * 1e6, 3.3675 * 1e6 };
		mu = { 0.016, 0.07215 };
		omega = { 0.01142, 0.2514 };
		rho_cr = { 162.66, 232 };
		s = { -0.154 , -0.04183 };
		k = { {0., 0.03}, { 0.03, 0.} };
		gamma = { 0.432, 1. };
	}
	case 3: // metan, pentan, etan
	{
		Tcr = { 190.564, 469.65, 305.51 };
		p_cr = { 4.5992 * 1e6, 3.3675 * 1e6, 4.8711 * 1e6 };
		mu = { 0.016, 0.07215, 0.03007 };
		omega = { 0.01142, 0.2514, 0.0979 };
		rho_cr = { 162.66, 232, 200 };
		s = { -0.154 , -0.04183, 0.1002 };
		k = { {0, 0.03, 0.005}, { 0.03, 0, 0.01}, {0.005, 0.01, 0} };
		gamma = { 0.432, 1., 0.5 };
	}
	};
	dx = { 0, 1, -1, 0,  0, 0,  0, 1, -1, -1,  1, 1, -1, -1,  1, 0,  0,  0,  0 };
	dy = { 0, 0,  0, 1, -1, 0,  0, 1,  1, -1, -1, 0,  0,  0,  0, 1, -1, -1,  1 };
	dz = { 0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1,  1, -1, -1, 1,  1, -1, -1 };
	index = { 0, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
	opposite_index = { 0, 2, 1, 4, 3, 6, 5, 9, 10, 7,  8, 13, 14, 11, 12, 17, 18, 15, 16 };
}

inline double scalar_product(vector<double> vec1, double vec2, double vec3, double vec4) {
	return vec1[0] * vec2 + vec1[1] * vec3 + vec1[2] * vec4;
}

inline double squaring(double vec1, double vec2, double vec3) {
	return vec1 * vec1 + vec2 * vec2 + vec3 * vec3;
}

double LatticeBoltzmannComponents::F_e(vector<double> vec1, double vec2, double vec3, double vec4, double w, double rho)
{
	return w * rho* (1 + scalar_product(vec1, vec2, vec3, vec4) / (teta)+scalar_product(vec1, vec2, vec3, vec4) * scalar_product(vec1, vec2, vec3, vec4) / (2. * teta * teta) - squaring(vec2, vec3, vec4) / 2. / (teta));
}

double LatticeBoltzmannComponents::F(double f, double f_eq, double f_eq1, double f_eq2)
{
	return f + (f_eq - f) / tau + (f_eq1 - f_eq2);
}

double LatticeBoltzmannComponents::a(const double temperature, double omega)
{
	double m = 0.382144 + 1.476905 * omega - 0.134488 * omega * omega;
	double a = pow((1 + m * (1 - sqrt(temperature))), 2);
	return a;
}

double LatticeBoltzmannComponents::PressurePengRobinson(double rho, const double temperature, double omega)
{
	double pressure = 1 / 0.307 * (temperature / (1. / rho - 0.253) -
		1.487 * a(temperature, omega) / (1. / rho / rho + 2 * 0.253 / rho - 0.253 * 0.253));
	return pressure;
}

double LatticeBoltzmannComponents::PressurePengRobinsonMultyComponent(vector<double> Rho)
{
	double B = 0, A = 0, D = 0, S = 0;
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

void LatticeBoltzmannComponents::SaveVTKFile(int tStep)
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

	vtk_file << "SCALARS rho_sum double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int l = 1; l < N_z + 1; l++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int i = 1; i < N_x + 1; i++) {
				double sum = 0;
				for (int numComp = 0; numComp < numberComponent; numComp++) {
					sum += rho[numComp][i][j][l];
				}
				vtk_file << sum << " ";
			}
		}
	}
	vtk_file << endl;

	for (int number = 1; number < numberComponent + 1; number++) {
		vtk_file << "VECTORS uflow" << number <<" double\n";
		for (int l = 1; l < N_z + 1; l++)
			for (int j = 1; j < N_y + 1; j++)
				for (int i = 1; i < N_x + 1; i++)
					vtk_file << ux[number - 1][i][j][l] + g / 2 << " " << uy[number - 1][i][j][l] << " " << uz[number - 1][i][j][l] << " ";
		vtk_file << endl;

		vtk_file << "SCALARS rho" << number << " double 1\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int l = 1; l < N_z + 1; l++)
			for (int j = 1; j < N_y + 1; j++)
				for (int i = 1; i < N_x + 1; i++)
					vtk_file << rho[number - 1][i][j][l] << " ";
		vtk_file << endl;
	}

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

void LatticeBoltzmannComponents::initialize()
{
	/*vector<double> temp(3, 0.0);
	for (int j = 1; j < N_y + 1; j++) {
		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				temp = mixture(percent1, percent2);
				rho[0][i][j][l] = temp[0];
				rho[1][i][j][l] = temp[1] + 2 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
				rho[2][i][j][l] = temp[2];
			}
		}
	}*/
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
	
	/* vector of possible velocities*/
	for (int j = 1; j < 19; j++) {
		c[j][0] = h * dx[j] / delta_t;
		c[j][1] = h * dy[j] / delta_t;
		c[j][2] = h * dz[j] / delta_t;
	}

	/*coefficients w_k */
	w[0] = 1. / 3.;

	for (int i = 1; i < 7; i++) {
		w[i] = 1. / 18.;
	}

	for (int i = 7; i < 19; i++) {
		w[i] = 1. / 36.;
	}

	/* coefficients G_k*/
	for (int i = 0; i < 7; i++) {
		G[i] = 1.;
	}

	for (int i = 7; i < 19; i++) {
		G[i] = 1. / 2.;
	}

#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				mask[i][j][l] = 0;
				for (int s = 0; s < 19; s++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						f[numComp][s][i][j][l] = F_e(c[s], ux[numComp][i][j][l], uy[numComp][i][j][l], uz[numComp][i][j][l], w[s], rho[numComp][i][j][l]);
					}
				}
			}
		}
	}
}

vector<double> LatticeBoltzmannComponents::mixture(double per1, double per2)
{
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

// мб что-то с мрт скопировала, если не будет работать, проверить
void LatticeBoltzmannComponents::setPeriodCondition( int numComp)
{
	X_condition[0] = true;
	X_condition[1] = false;
	X_condition[2] = false;
	wall_condition[0] = true;
	wall_condition[1] = false;
	wall_condition[2] = false;
	for (size_t l = 1; l < N_z + 1; l++) {
		for (size_t j = 1; j < N_y + 1; j++) {
			buf[numComp][1][0][j][l] = buf[numComp][1][N_x][j][l];
			buf[numComp][2][N_x + 1][j][l] = buf[numComp][2][1][j][l];
		};
	};
	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 1; l < N_z + 1; l++) {
			buf[numComp][7][0][j][l] = buf[numComp][7][N_x][j + 1][l];
		}
	}
	for (size_t j = 1; j < N_y; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][8][N_x + 1][j][l] = buf[numComp][8][1][j + 1][l];
		}
	}
	for (size_t j = 2; j < N_y + 1; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][9][N_x + 1][j][l] = buf[numComp][9][1][j - 1][l];
		}
	}
	for (size_t j = 2; j < N_y + 1; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][10][0][j][l] = buf[numComp][10][N_x][j - 1][l];
		}
	}

	for (size_t l = 1; l < N_z; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][11][0][j][l] = buf[numComp][11][N_x][j][l + 1];
		}
	}

	for (size_t l = 1; l < N_z; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][12][N_x + 1][j][l] = buf[numComp][12][1][j][l + 1];
		}
	}

	for (size_t l = 2; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][13][N_x + 1][j][l] = buf[numComp][13][1][j][l - 1];
		}
	}

	for (size_t l = 2; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][14][0][j][l] = buf[numComp][14][N_x][j][l - 1];
		}
	}
		
	// period walls
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			buf[numComp][5][i][j][0] = buf[numComp][5][i][j][N_z];
			buf[numComp][6][i][j][N_z + 1] = buf[numComp][6][i][j][1];
			buf[numComp][11][i][j][0] = buf[numComp][11][i][j][N_z];
			buf[numComp][12][i][j][0] = buf[numComp][12][i][j][N_z];
			buf[numComp][13][i][j][N_z + 1] = buf[numComp][13][i][j][1];
			buf[numComp][14][i][j][N_z + 1] = buf[numComp][14][i][j][1];
			buf[numComp][15][i][j][0] = buf[numComp][15][i][j][N_z];
			buf[numComp][16][i][j][0] = buf[numComp][16][i][j][N_z];
			buf[numComp][17][i][j][N_z + 1] = buf[numComp][17][i][j][1];
			buf[numComp][18][i][j][N_z + 1] = buf[numComp][18][i][j][1];
		}
	}

	for (int i = 1; i < N_x + 1; i++) {
		for (int l = 1; l < N_z + 1; l++) {
			buf[numComp][3][i][0][l] = buf[numComp][3][i][N_y][l];
			buf[numComp][4][i][N_y + 1][l] = buf[numComp][4][i][1][l];
			buf[numComp][7][i][0][l] = buf[numComp][7][i][N_y][l];
			buf[numComp][8][i][0][l] = buf[numComp][8][i][N_y][l];
			buf[numComp][9][i][N_y + 1][l] = buf[numComp][9][i][1][l];
			buf[numComp][10][i][N_y + 1][l] = buf[numComp][10][i][1][l];
			buf[numComp][15][i][0][l] = buf[numComp][15][i][N_y][l];
			buf[numComp][16][i][N_y + 1][l] = buf[numComp][16][i][1][l];
			buf[numComp][17][i][N_y + 1][l] = buf[numComp][17][i][1][l];
			buf[numComp][18][i][0][l] = buf[numComp][18][i][N_y][l];
		}
	}

	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 1; l < N_z + 1; l++) {
			buf[numComp][1][0][j][l] = buf[numComp][1][N_x][j][l];
			buf[numComp][2][N_x + 1][j][l] = buf[numComp][2][1][j][l];
			buf[numComp][7][0][j][l] = buf[numComp][7][N_x][j][l];
			buf[numComp][8][N_x + 1][j][l] = buf[numComp][8][1][j][l];
			buf[numComp][9][N_x + 1][j][l] = buf[numComp][9][1][j][l];
			buf[numComp][10][0][j][l] = buf[numComp][10][N_x][j][l];
			buf[numComp][11][0][j][l] = buf[numComp][11][N_x][j][l];
			buf[numComp][12][N_x + 1][j][l] = buf[numComp][12][1][j][l];
			buf[numComp][13][N_x + 1][j][l] = buf[numComp][13][1][j][l];
			buf[numComp][14][0][j][l] = buf[numComp][14][N_x][j][l];
		}
	}

	for (int i = 1; i < N_x + 1; i++) {
		buf[numComp][15][i][0][0] = buf[numComp][15][i][N_y][N_z];
		buf[numComp][17][i][N_y + 1][N_z + 1] = buf[numComp][17][i][1][1];
		buf[numComp][16][i][N_y + 1][0] = buf[numComp][16][i][1][N_z];
		buf[numComp][18][i][0][N_z + 1] = buf[numComp][18][i][N_y][1];
	}

	for (int j = 1; j < N_y + 1; j++) {
		buf[numComp][11][0][j][0] = buf[numComp][11][N_x][j][N_z];
		buf[numComp][12][N_x + 1][j][0] = buf[numComp][12][1][j][N_z];
		buf[numComp][13][N_x + 1][j][N_z + 1] = buf[numComp][13][1][j][1];
		buf[numComp][14][0][j][N_z + 1] = buf[numComp][14][N_x][j][1];
	}

	for (int l = 1; l < N_z + 1; l++) {
		buf[numComp][7][0][0][l] = buf[numComp][7][N_x][N_y][l];
		buf[numComp][8][N_x + 1][0][l] = buf[numComp][8][1][N_y][l];
		buf[numComp][9][N_x + 1][N_y + 1][l] = buf[numComp][9][1][1][l];
		buf[numComp][10][0][N_y + 1][l] = buf[numComp][10][N_x][1][l];
	}
	
}

void LatticeBoltzmannComponents::setFreeConditionByX( int numComp)
{
	X_condition[0] = false;
	X_condition[1] = true;
	X_condition[2] = false;
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

void LatticeBoltzmannComponents::setWalls( int numComp)
{
	wall_condition[0] = false;
	wall_condition[1] = false;
	wall_condition[2] = true;
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

void LatticeBoltzmannComponents::setWallsByX( int numComp)
{
	X_condition[0] = false;
	X_condition[1] = false;
	X_condition[2] = true;

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

void LatticeBoltzmannComponents::setFiCondition()
{
	if (X_condition[0] == true && wall_condition[0] == true) 
	{
		// period condition
		Fi[0][0][0] = wetWalls * Fi[N_x][N_y][N_z];
		Fi[0][N_y + 1][0] = wetWalls * Fi[N_x][1][N_z];
		Fi[0][0][N_z + 1] = wetWalls * Fi[N_x][N_y][1];
		Fi[0][N_y + 1][N_z + 1] = wetWalls * Fi[N_x][1][1];
		Fi[N_x + 1][0][0] = wetWalls * Fi[1][N_y][N_z];
		Fi[N_x + 1][N_y + 1][0] = wetWalls * Fi[1][1][N_z];
		Fi[N_x + 1][0][N_z + 1] = wetWalls * Fi[1][N_y][1];
		Fi[N_x + 1][N_y + 1][N_z + 1] = wetWalls * Fi[1][1][1];

		for (int l = 1; l < N_z + 1; l++) {
			Fi[0][0][l] = wetWalls * Fi[N_x][N_y][l];
			Fi[0][N_y + 1][l] = wetWalls * Fi[N_x][1][l];
			Fi[N_x + 1][0][l] = wetWalls * Fi[1][N_y][l];
			Fi[N_x + 1][N_y + 1][l] = wetWalls * Fi[1][1][l];
		}

		for (int i = 1; i < N_x + 1; i++) {
			Fi[i][0][0] = wetWalls * Fi[i][N_y][N_z];
			Fi[i][N_y + 1][0] = wetWalls * Fi[i][1][N_z];
			Fi[i][0][N_z + 1] = wetWalls * Fi[i][N_y][1];
			Fi[i][N_y + 1][N_z + 1] = wetWalls * Fi[i][1][1];
		}

		for (int j = 1; j < N_y + 1; j++) {
			Fi[0][j][0] = wetWalls * Fi[N_x][j][N_z];
			Fi[N_x + 1][j][0] = wetWalls * Fi[1][j][N_z];
			Fi[0][j][N_z + 1] = wetWalls * Fi[N_x][j][1];
			Fi[N_x + 1][j][N_z + 1] = wetWalls * Fi[1][j][1];
		}

		for (int i = 1; i < N_x + 1; i++) {
			for (int j = 1; j < N_y + 1; j++) {
				Fi[i][j][0] = wetWalls * Fi[i][j][N_z];
				Fi[i][j][N_z + 1] = wetWalls * Fi[i][j][1];
			}
		}

		for (int i = 1; i < N_x + 1; i++) {
			for (int l = 1; l < N_z + 1; l++) {
				Fi[i][0][l] = wetWalls * Fi[i][N_y][l];
				Fi[i][N_y + 1][l] = wetWalls * Fi[i][1][l];
			}
		}

		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				Fi[0][j][l] = wetWalls * Fi[N_x][j][l];
				Fi[N_x + 1][j][l] = wetWalls * Fi[1][j][l];
			}
		}
	}

	if(wall_condition[2] == true)
	{
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
	}

	if (X_condition[2] == true)
	{
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				Fi[0][j][l] = wetWalls * Fi[1][j][l];
				Fi[N_x + 1][j][l] = wetWalls * Fi[N_x][j][l];
			}
		}
	}

	if (X_condition[1] == true)
	{
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				Fi[0][j][l] = Fi[1][j][l];
				Fi[N_x + 1][j][l] = Fi[N_x][j][l];
			}
		}
	}
}


void LatticeBoltzmannComponents::TimeStep(int tStep)
{
	buf = f;

  if (tStep < RelaxTime) {
	  for (int numComp = 0; numComp < numberComponent; numComp++) {
		  setWalls(numComp);
		  setWallsByX(numComp);
		  //setPeriodCondition(buf);
	  }
	}
	else {
	  g = 1e-1;
	  for (int numComp = 0; numComp < numberComponent; numComp++) {
		  setFreeConditionByX(numComp);
		  setWalls(numComp);
	  }
	}

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
	for (int number = 1; number < numberComponent + 1; number++){
		Full_rho[number - 1] = 0.;
		rho_min[number - 1] = 1000.;
		rho_max[number - 1] = -1000.;
	}
	rho_mix_max = -1000;
	rho_mix_min = 1000;

//#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				for (int number = 1; number < numberComponent + 1; number++) {
					Full_rho[number - 1] += rho[number - 1][i][j][l];
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
				for (int number = 1; number < numberComponent + 1; number++) {
					if (rho_min[number - 1] > rho[number - 1][i][j][l]) rho_min[number - 1] = rho[number - 1][i][j][l];
					if (rho_max[number - 1] < rho[number - 1][i][j][l]) rho_max[number - 1] = rho[number - 1][i][j][l];
				}
				if (rho_mix_min > Sum_rho[i][j][l]) rho_mix_min = Sum_rho[i][j][l];
				if (rho_mix_max < Sum_rho[i][j][l]) rho_mix_max = Sum_rho[i][j][l];

			}
		};
	};

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
						ux[numComp][i][j][l] += f[numComp][s][i][j][l] * c[s][0] / rho[numComp][i][j][l];
						uy[numComp][i][j][l] += f[numComp][s][i][j][l] * c[s][1] / rho[numComp][i][j][l];
						uz[numComp][i][j][l] += f[numComp][s][i][j][l] * c[s][2] / rho[numComp][i][j][l];
					}
				}
			}
		}
	}

	/* pressure */
	vector<double> temp(numberComponent, 0.0);
#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				for (int numComp = 0; numComp < numberComponent; numComp++) {
					temp[numComp] = rho[numComp][i][j][l];
				}
				pressure[i][j][l] = PressurePengRobinsonMultyComponent(temp);
			}
		}
	}

	/* "effective" density */
#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				Fi[i][j][l] = sqrt(Sum_rho[i][j][l] * teta - pressure[i][j][l]);
			}
		}
	}

	setFiCondition();

#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				if (mask[i][j][l] == 1.0) {
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
							dux[numComp][i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dx[s] +
								A * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dx[s]) * delta_t * gamma[numComp] / gamma_rho[i][j][l] / h;
							duy[numComp][i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dy[s] +
								A * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dy[s]) * delta_t * gamma[numComp] / gamma_rho[i][j][l] / h;
							duz[numComp][i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dz[s] +
								A * G[s] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * Fi[i + dx[s]][j + dy[s]][l + dz[s]] * dz[s]) * delta_t * gamma[numComp] / gamma_rho[i][j][l] / h;
						}
					}
				}
			}
		}
	}

	/* the law of conservation of impulse */
	for (int number = 1; number < numberComponent + 1; number++) {
		Full_velocity_x[number - 1] = 0;
		Full_velocity_y[number - 1] = 0;
		Full_velocity_z[number - 1] = 0;
	}
//#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				for (int numComp = 0; numComp < numberComponent; numComp++) {
					if (mask[i][j][l] == 0) {
						Full_velocity_x[numComp] += rho[numComp][i][j][l] * (ux[numComp][i][j][l] + (dux[numComp][i][j][l] + g) / 2);
						Full_velocity_y[numComp] += rho[numComp][i][j][l] * (uy[numComp][i][j][l] + duy[numComp][i][j][l] * 0.5);
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
}

void LatticeBoltzmannComponents::collisions()
{
#pragma omp parallel for
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				for (int s = 0; s < 19; s++) {
					for (int numComp = 0; numComp < numberComponent; numComp++) {
						f[numComp][s][i][j][l] = F(f[numComp][s][i][j][l],
							F_e(c[s], ux_all[i][j][l], uy_all[i][j][l], uz_all[i][j][l], w[s], rho[numComp][i][j][l]),
							F_e(c[s], ux[numComp][i][j][l] + dux[numComp][i][j][l] + g, uy[numComp][i][j][l] + duy[numComp][i][j][l], uz[numComp][i][j][l] + duz[numComp][i][j][l], w[s], rho[numComp][i][j][l]),
							F_e(c[s], ux[numComp][i][j][l], uy[numComp][i][j][l], uz[numComp][i][j][l], w[s], rho[numComp][i][j][l]));
					}
				}
			};
		};
	};
}

void LatticeBoltzmannComponents::data(int tStep)
{
	if (tStep % 100 == 0)
	{
		SaveVTKFile(tStep);
		cout << "Summa mass 1 = " << Full_rho[0] << endl;
		if (numberComponent > 2) {
			cout << "Summa mass 2 = " << Full_rho[1] << endl;
		}
		if (numberComponent == 3) {
			cout << "Summa mass 3 = " << Full_rho[2] << endl;
		}
		if (numberComponent == 1) {
			cout << " Impulse for " << tStep << " step (x) = " << Full_velocity_x[0] << endl;
			cout << " Impulse for " << tStep << " step (y) = " << Full_velocity_y[0] << endl;
			cout << " Impulse for " << tStep << " step (z) = " << Full_velocity_z[0] << endl;
		}
		if (numberComponent == 2) {
			cout << " Impulse for " << tStep << " step (x) = " << Full_velocity_x[0] << " and  " << Full_velocity_x[1] << endl;
			cout << " Impulse for " << tStep << " step (y) = " << Full_velocity_y[0] << " and  " << Full_velocity_y[1] << endl;
			cout << " Impulse for " << tStep << " step (z) = " << Full_velocity_z[0] << " and  " << Full_velocity_z[1] << endl;
		}
		if (numberComponent == 3) {
			cout << " Impulse for " << tStep << " step (x) = " << Full_velocity_x[0] << " and  " << Full_velocity_x[1] << " and  " << Full_velocity_x[2] << endl;
			cout << " Impulse for " << tStep << " step (y) = " << Full_velocity_y[0] << " and  " << Full_velocity_y[1] << " and  " << Full_velocity_y[2] << endl;
			cout << " Impulse for " << tStep << " step (z) = " << Full_velocity_z[0] << " and  " << Full_velocity_z[1] << " and  " << Full_velocity_z[2] << endl;
		}

		cout << " rho mix min = " << rho_mix_min << " and rho mix max = " << rho_mix_max << endl;
		cout << " rho1 min = " << rho_min[0] << " and rho1 max = " << rho_max[0] << endl;
		if (numberComponent >= 2) {
			cout << " rho2 min = " << rho_min[1] << " and rho2 max = " << rho_max[1] << endl;
		}
		if (numberComponent >= 3) {
			cout << " rho3 min = " << rho_min[2] << " and rho3 max = " << rho_max[2] << endl;
		}
		cout << g << endl;
	}
}

MRT::MRT():
m(vector<double>(kMax)),
m_eq(vector<double>(kMax))
{
	dx = { 0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0,  0 };
	dy = { 0, 0,  0, 1, -1, 0,  0, 1,  1, -1, -1, 0,  0,  0,  0, 1, -1,  1, -1 };
	dz = { 0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1,  1, -1, -1, 1,  1, -1, -1 };
	index = {          0, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
	opposite_index = { 0, 2, 1, 4, 3, 6, 5, 10, 9, 8,  7, 14, 13, 12, 11, 18, 17, 16, 15 };
	s_ii = { 0, s1, s2, 0,  s4, 0,  s4, 0, s4,  s9, s10, s9, s10, s13, s13, s13, s16, s16, s16 };
}

void MRT::setPeriodCondition(int numComp)
{
	X_condition[0] = true;
	X_condition[1] = false;
	X_condition[2] = false;
	wall_condition[0] = true;
	wall_condition[1] = false;
	wall_condition[2] = false;
	for (size_t l = 1; l < N_z + 1; l++) {
		for (size_t j = 1; j < N_y + 1; j++) {
			buf[numComp][1][0][j][l] = buf[numComp][1][N_x][j][l];
			buf[numComp][2][N_x + 1][j][l] = buf[numComp][2][1][j][l];
		};
	};
	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 1; l < N_z + 1; l++) {
			buf[numComp][7][0][j][l] = buf[numComp][7][N_x][j + 1][l];
		}
	}
	for (size_t j = 1; j < N_y; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][8][N_x + 1][j][l] = buf[numComp][8][1][j + 1][l];
		}
	}
	for (size_t j = 2; j < N_y + 1; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][10][N_x + 1][j][l] = buf[numComp][10][1][j - 1][l];
		}
	}
	for (size_t j = 2; j < N_y + 1; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][9][0][j][l] = buf[numComp][9][N_x][j - 1][l];
		}
	}

	for (size_t l = 1; l < N_z; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][11][0][j][l] = buf[numComp][11][N_x][j][l + 1];
		}
	}

	for (size_t l = 1; l < N_z; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][12][N_x + 1][j][l] = buf[numComp][12][1][j][l + 1];
		}
	}

	for (size_t l = 2; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][14][N_x + 1][j][l] = buf[numComp][14][1][j][l - 1];
		}
	}

	for (size_t l = 2; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][13][0][j][l] = buf[numComp][13][N_x][j][l - 1];
		}
	}

	// period walls
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			buf[numComp][5][i][j][0] = buf[numComp][5][i][j][N_z];
			buf[numComp][6][i][j][N_z + 1] = buf[numComp][6][i][j][1];
			buf[numComp][11][i][j][0] = buf[numComp][11][i][j][N_z];
			buf[numComp][12][i][j][0] = buf[numComp][12][i][j][N_z];
			buf[numComp][14][i][j][N_z + 1] = buf[numComp][14][i][j][1];
			buf[numComp][13][i][j][N_z + 1] = buf[numComp][13][i][j][1];
			buf[numComp][15][i][j][0] = buf[numComp][15][i][j][N_z];
			buf[numComp][16][i][j][0] = buf[numComp][16][i][j][N_z];
			buf[numComp][18][i][j][N_z + 1] = buf[numComp][18][i][j][1];
			buf[numComp][17][i][j][N_z + 1] = buf[numComp][17][i][j][1];
		}
	}

	for (int i = 1; i < N_x + 1; i++) {
		for (int l = 1; l < N_z + 1; l++) {
			buf[numComp][3][i][0][l] = buf[numComp][3][i][N_y][l];
			buf[numComp][4][i][N_y + 1][l] = buf[numComp][4][i][1][l];
			buf[numComp][7][i][0][l] = buf[numComp][7][i][N_y][l];
			buf[numComp][8][i][0][l] = buf[numComp][8][i][N_y][l];
			buf[numComp][10][i][N_y + 1][l] = buf[numComp][10][i][1][l];
			buf[numComp][9][i][N_y + 1][l] = buf[numComp][9][i][1][l];
			buf[numComp][15][i][0][l] = buf[numComp][15][i][N_y][l];
			buf[numComp][16][i][N_y + 1][l] = buf[numComp][16][i][1][l];
			buf[numComp][18][i][N_y + 1][l] = buf[numComp][18][i][1][l];
			buf[numComp][17][i][0][l] = buf[numComp][17][i][N_y][l];
		}
	}

	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 1; l < N_z + 1; l++) {
			buf[numComp][1][0][j][l] = buf[numComp][1][N_x][j][l];
			buf[numComp][2][N_x + 1][j][l] = buf[numComp][2][1][j][l];
			buf[numComp][7][0][j][l] = buf[numComp][7][N_x][j][l];
			buf[numComp][8][N_x + 1][j][l] = buf[numComp][8][1][j][l];
			buf[numComp][10][N_x + 1][j][l] = buf[numComp][10][1][j][l];
			buf[numComp][9][0][j][l] = buf[numComp][9][N_x][j][l];
			buf[numComp][11][0][j][l] = buf[numComp][11][N_x][j][l];
			buf[numComp][12][N_x + 1][j][l] = buf[numComp][12][1][j][l];
			buf[numComp][14][N_x + 1][j][l] = buf[numComp][14][1][j][l];
			buf[numComp][13][0][j][l] = buf[numComp][13][N_x][j][l];
		}
	}

	for (int i = 1; i < N_x + 1; i++) {
		buf[numComp][15][i][0][0] = buf[numComp][15][i][N_y][N_z];
		buf[numComp][18][i][N_y + 1][N_z + 1] = buf[numComp][18][i][1][1];
		buf[numComp][16][i][N_y + 1][0] = buf[numComp][16][i][1][N_z];
		buf[numComp][17][i][0][N_z + 1] = buf[numComp][17][i][N_y][1];
	}

	for (int j = 1; j < N_y + 1; j++) {
		buf[numComp][11][0][j][0] = buf[numComp][11][N_x][j][N_z];
		buf[numComp][12][N_x + 1][j][0] = buf[numComp][12][1][j][N_z];
		buf[numComp][14][N_x + 1][j][N_z + 1] = buf[numComp][14][1][j][1];
		buf[numComp][13][0][j][N_z + 1] = buf[numComp][13][N_x][j][1];
	}

	for (int l = 1; l < N_z + 1; l++) {
		buf[numComp][7][0][0][l] = buf[numComp][7][N_x][N_y][l];
		buf[numComp][8][N_x + 1][0][l] = buf[numComp][8][1][N_y][l];
		buf[numComp][10][N_x + 1][N_y + 1][l] = buf[numComp][10][1][1][l];
		buf[numComp][9][0][N_y + 1][l] = buf[numComp][9][N_x][1][l];
	}
}

void MRT::setFreeConditionByX(int numComp)
{
	X_condition[0] = false;
	X_condition[1] = true;
	X_condition[2] = false;
	buf[numComp][1][0] = buf[numComp][1][1];
	buf[numComp][2][N_x + 1] = buf[numComp][2][N_x];
	buf[numComp][7][0] = buf[numComp][7][1];
	buf[numComp][8][N_x + 1] = buf[numComp][8][N_x];
	buf[numComp][10][N_x + 1] = buf[numComp][10][N_x];
	buf[numComp][9][0] = buf[numComp][9][1];
	buf[numComp][11][0] = buf[numComp][11][1];
	buf[numComp][12][N_x + 1] = buf[numComp][12][N_x];
	buf[numComp][14][N_x + 1] = buf[numComp][14][N_x];
	buf[numComp][13][0] = buf[numComp][13][1];
}

void MRT::setWalls(int numComp)
{
	wall_condition[0] = false;
	wall_condition[1] = false;
	wall_condition[2] = true;
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
			buf[numComp][7][i][0][l] = buf[numComp][10][i + 1][1][l];
			buf[numComp][9][i][N_y + 1][l] = buf[numComp][8][i + 1][N_y][l];
		}
	}
	for (int i = 1; i < N_x + 2; i++) {
		for (int l = 1; l < N_z + 1; l++) {
			buf[numComp][10][i][N_y + 1][l] = buf[numComp][7][i - 1][N_y][l];
			buf[numComp][8][i][0][l] = buf[numComp][9][i - 1][1][l];
		}
	}

	for (int i = 0; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			buf[numComp][11][i][j][0] = buf[numComp][14][i + 1][j][1];
			buf[numComp][13][i][j][N_z + 1] = buf[numComp][12][i + 1][j][N_z];
		}
	}
	for (int i = 1; i < N_x + 2; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			buf[numComp][14][i][j][N_z + 1] = buf[numComp][11][i - 1][j][N_z];
			buf[numComp][12][i][j][0] = buf[numComp][13][i - 1][j][1];
		}
	}

	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 0; j < N_y + 1; j++) {
			buf[numComp][15][i][j][0] = buf[numComp][18][i][j + 1][1];
			buf[numComp][17][i][j][N_z + 1] = buf[numComp][16][i][j + 1][N_z];
		}
	}
	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 2; j++) {
			buf[numComp][18][i][j][N_z + 1] = buf[numComp][15][i][j - 1][N_z];
			buf[numComp][16][i][j][0] = buf[numComp][17][i][j - 1][1];
		}
	}

	for (int i = 1; i < N_x + 1; i++) {
		for (int l = 0; l < N_z + 1; l++) {
			buf[numComp][15][i][0][l] = buf[numComp][18][i][1][l + 1];
			buf[numComp][16][i][N_y + 1][l] = buf[numComp][17][i][N_y][l + 1];
		}
	}
	for (int i = 0; i < N_x + 1; i++) {
		for (int l = 1; l < N_z + 2; l++) {
			buf[numComp][18][i][N_y + 1][l] = buf[numComp][15][i][N_y][l - 1];
			buf[numComp][17][i][0][l] = buf[numComp][16][i][1][l - 1];
		}
	}
}

void MRT::setWallsByX(int numComp)
{
	X_condition[0] = false;
	X_condition[1] = false;
	X_condition[2] = true;

	for (size_t l = 1; l < N_z + 1; l++) {
		for (size_t j = 1; j < N_y + 1; j++) {
			buf[numComp][1][0][j][l] = buf[numComp][2][1][j][l];
			buf[numComp][2][N_x + 1][j][l] = buf[numComp][1][N_x][j][l];
		};
	};
	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 1; l < N_z + 1; l++) {
			buf[numComp][7][0][j][l] = buf[numComp][10][1][j + 1][l];
		}
	}
	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][8][N_x + 1][j][l] = buf[numComp][9][N_x][j + 1][l];
		}
	}
	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][10][N_x + 1][j][l] = buf[numComp][7][N_x][j - 1][l];
		}
	}
	for (size_t j = 1; j < N_y + 1; j++) {
		for (size_t l = 0; l < N_z + 1; l++) {
			buf[numComp][9][0][j][l] = buf[numComp][8][1][j - 1][l];
		}
	}

	for (size_t l = 1; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][11][0][j][l] = buf[numComp][14][1][j][l + 1];
		}
	}

	for (size_t l = 1; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][12][N_x + 1][j][l] = buf[numComp][13][N_x][j][l + 1];
		}
	}

	for (size_t l = 1; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][14][N_x + 1][j][l] = buf[numComp][11][N_x][j][l - 1];
		}
	}

	for (size_t l = 1; l < N_z + 1; l++) {
		for (size_t j = 0; j < N_y + 1; j++) {
			buf[numComp][13][0][j][l] = buf[numComp][12][1][j][l - 1];
		}
	}
}

double MRT::Kroneker(int a, int b)
{
	double number;
	if (a == b) {
		number = 1.;
	}
	else {
		number = 0.;
	}
	return number;
}

void MRT::collisions()
{
	vector<double> moment(19);
#pragma omp parallel for

	for (int i = 1; i < N_x + 1; i++) {
		for (int j = 1; j < N_y + 1; j++) {
			for (int l = 1; l < N_z + 1; l++) {
				for (int numComp = 0; numComp < numberComponent; numComp++) {
					for (int s = 0; s < kMax; s++) {
						m[s] = 0;
						for (int s1 = 0; s1 < kMax; s1++) {
							m[s] += M[s][s1] * f[numComp][s1][i][j][l];
						}
					}
					m_eq[0] = rho[numComp][i][j][l];
					m_eq[1] = -11 * rho[numComp][i][j][l] + 19 * rho[numComp][i][j][l] * squaring(ux_all[i][j][l], uy_all[i][j][l], uz_all[i][j][l]) * delta_t * delta_t / h / h;
					m_eq[2] = w_e * rho[numComp][i][j][l] + w_ej * rho[numComp][i][j][l] * squaring(ux_all[i][j][l], uy_all[i][j][l], uz_all[i][j][l]) * delta_t * delta_t / h / h;
					m_eq[3] = ux_all[i][j][l] * rho[numComp][i][j][l] * delta_t / h;
					m_eq[4] = -2. / 3. * ux_all[i][j][l] * rho[numComp][i][j][l] * delta_t / h;
					m_eq[5] = uy_all[i][j][l] * rho[numComp][i][j][l] * delta_t / h;
					m_eq[6] = -2. / 3. * uy_all[i][j][l] * rho[numComp][i][j][l] * delta_t / h;
					m_eq[7] = uz_all[i][j][l] * rho[numComp][i][j][l] * delta_t / h;
					m_eq[8] = -2. / 3. * uz_all[i][j][l] * rho[numComp][i][j][l] * delta_t / h;
					m_eq[9] = 1 / rho[numComp][i][j][l] * (2 * m_eq[3] * m_eq[3] - (m_eq[5] * m_eq[5] + m_eq[7] * m_eq[7]));
					m_eq[10] = w_xx * m_eq[9];
					m_eq[11] = 1 / rho[numComp][i][j][l] * (m_eq[5] * m_eq[5] - m_eq[7] * m_eq[7]);
					m_eq[12] = w_xx * m_eq[11];
					m_eq[13] = 1 / rho[numComp][i][j][l] * m_eq[3] * m_eq[5];
					m_eq[14] = 1 / rho[numComp][i][j][l] * m_eq[7] * m_eq[5];
					m_eq[15] = 1 / rho[numComp][i][j][l] * m_eq[3] * m_eq[7];
					m_eq[16] = 0;
					m_eq[17] = 0;
					m_eq[18] = 0;

					for (int s = 0; s < kMax; s++) {
						moment[s] = 0;
						for (int s1 = 0; s1 < kMax; s1++) {
							for (int s2 = 0; s2 < kMax; s2++) {
								moment[s] += -Mi[s][s1] * s_ii[s2] * (m[s1] - m_eq[s1]) * Kroneker(s1, s2);
							}
						}
					}
					for (int s = 0; s < kMax; s++) {
						f[numComp][s][i][j][l] = f[numComp][s][i][j][l] + moment[s] + F_e(c[s], ux[numComp][i][j][l] + dux[numComp][i][j][l] + g, uy[numComp][i][j][l] + duy[numComp][i][j][l],
							uz[numComp][i][j][l] + duz[numComp][i][j][l], w[s], rho[numComp][i][j][l]) -
							F_e(c[s], ux[numComp][i][j][l], uy[numComp][i][j][l], uz[numComp][i][j][l], w[s], rho[numComp][i][j][l]);
					}
				}
			}
		};
	};
}

