#pragma once
#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
#include <omp.h>

using namespace std;

class LatticeBoltzmannComponents {
public:
	LatticeBoltzmannComponents();

	/* равновесные функции распределения */
	double F_e(vector<double> vec1, double vec2, double vec3, double vec4, double w, double rho);

	/* решеточное уравнение Больцмана */
	double F(double f, double f_eq, double f_eq1, double f_eq2);

	/* уравнение состояния Пенга-Робинсона */
	double a(const double temperature, double omega);
	double PressurePengRobinson(double rho, const double temperature, double omega);
	double PressurePengRobinsonMultyComponent(vector<double> Rho);

	void SaveVTKFile(int tStep);

	void initialize();
	vector<double> mixture(double per1, double per2);

	void setPeriodCondition( int numComp);
	void setFreeConditionByX( int numComp);
	void setWalls( int numComp);
	void setWallsByX( int numComp);
	void setFiCondition();

	void TimeStep(int tStep);

private:
	const double R = 8.31446;
	const int kMax = 19;
	double delta_t = 1e-9;
	double h = 1e-6;
	int N_x = 30;
	int N_y = 5;
	int N_z = 5;
	int FullTime = 10000;
	int RelaxTime = 300;
	double teta = 1. / 3. * h * h / delta_t / delta_t;
	double tau = 1.;
	double A = -0.58;
	double T = 370;
	double g = 0;
	int numberComponent = 3;
	vector<double> Full_velocity_x, Full_velocity_y, Full_velocity_z, Full_rho, omega, Tcr, p_cr, mu, rho_cr, s, gamma, rho_max, G, w, rho_min;
	double wetWalls = 1;
	double wetObst = 1.;
	double sumG = 0, rho_mix_max, rho_mix_min;
	int volumeObs = 0;
	double a_0 = 0.4572793, b_0 = 0.07780669, percent1 = 0.2, percent2 = 0.6, rho_mix = 280;
	vector<vector<double>> k, c;
	vector<vector<vector<vector<vector<double>>>>> f, buf;
	vector<vector<vector<vector<double>>>> ux, uy, uz, dux, duy, duz, rho;
	vector<vector<vector<double>>> Fi, pressure, gamma_rho, Sum_rho, ux_all, uy_all, uz_all;
	vector<vector<vector<int>>> mask;
	vector <int> dx, dy, dz, opposite_index, index;
	vector<bool> wall_condition, X_condition;
};

/*можно еще сюда 2д подкрутить, вообще код супер будет; ну и мрт само собой
еще надо доработать, чтобы можно было безболезненно менять число компонент, и подумать, как сделать в прицнипе с любым чилом компонент так

*/