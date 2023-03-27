#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <fstream>
#include <math.h>
#include <complex>
#include <mkl.h>
#include <omp.h>
#include "linalg.h"
#include "TwoCubitsKernel.h"

using namespace std;

using TYPE = double;

const TYPE PI = acos(-1.0);

int main() {
	double time = omp_get_wtime();
	TYPE Theta = 0.0018; // Измеренный угол поворота
	TYPE w0q1 = 5.12 * (2 * PI) * 1e9; // Собственная частота первого кубита
	TYPE mu1 = -0.353 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	TYPE w0q2 = 5.350 * (2 * PI) * 1e9; // Собственная частота второго кубита
	TYPE mu2 = -0.35 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	TYPE g = 0.02 * (2 * PI) * 1e9; // параметр взаимодействия между кубитами
	TYPE wt = 5.13 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	TYPE wc = 5.21 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	TYPE Tmax = 50 * 1e-9; // Время интегрирования
	TYPE dt = 1e-13;
	int M = 10000; // Максимальное число периодов внешнего управляющего поля
	TYPE tau = 4 * 1e-12; // Длительность импульса

	TwoCubitsKernel<TYPE> kernel(Theta, w0q1, mu1, w0q2, mu2, g, wt, tau, Tmax, dt, M);
	cout << "Q1 = control qubit, Q2 = target qubit" << '\n';
	cout << "freq (Q1) = " << wc << '\n';
	cout << "freq (Q2) = " << wt << '\n';
	cout << "A (Q1) = " << kernel.Acon << '\n';
	cout << "A (Q2) = " << kernel.Atar << '\n';
	cout << "theta = " << Theta << '\n';
	cout << "tstep = " << dt << '\n';
	cout << "length = " << floor(kernel.Tmax * 1e9) << " ns" << '\n';

	const int PRECISION = 10;
	cout.precision(PRECISION);
	ofstream _01, _02, _03, _04;
	_01.open("imp.txt");
	_01.precision(PRECISION);
	_02.open("imp2.txt");
	_02.precision(PRECISION);
	int phi = 0;
	cout << "A1_stride | A2_stride | W1 | W2 | W3 | W4 | W5 | W6 | W7 | W8 | W9\n";
	for (TYPE A1_stride = 0; A1_stride <= 0.25 * kernel.Acon; A1_stride += 0.25 * kernel.Acon) {
		for (TYPE A2_stride = 0; A2_stride <= 0.25 * kernel.Acon; A2_stride += 0.25 * kernel.Atar) {
			auto initial_state = kernel.getEigVector(0);
			auto final_state = kernel.calculateDynamic(initial_state, A1_stride, A2_stride, phi, false, _01, _02);
			auto pops = kernel.getPopulations(final_state);
			cout << A1_stride << ' ' << A2_stride << ' ';
			for (auto j : pops) {
				cout << j << ' ';
			}
			cout << '\n';
		}
	}
	
	time = omp_get_wtime() - time;
	cout << "TIME ELAPSED = " << time << '\n';
	return 0;
}



