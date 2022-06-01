#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <complex>
#include <algorithm>
#include <omp.h>
#include <map>
#include "GeneticAlgorithm.h"
#include <functional>

using namespace std;

const int TYPE = 3;

const double tstep = 5e-14;
const double w01 = 4 * 2 * PI * 1e9;
const double w12 = w01 - 0.25 * 2 * PI * 1e9;

const double wt = 25 * 2 * PI * 1e9;
const double w = 4e-12;
const double T = PI / wt * 2;
double Theta = 0.001;
const int NumberOfCycles = 1;
const double NeededAngle = 0.024;

Kernel<TYPE> kernel(tstep, w01, w12, wt, w, T);

void bruteForce() {
	const int LEN = 8;
	vector<int> seq(LEN);
	pair<double, double> res = { INT_MAX, INT_MAX };
	vector<int> bestSeq;
	int bestLen;
	double bestTheta, bestF;
	function<void(int)> calc = [&](int pos) {
		if (pos == LEN) {
			vector<int> to_calc;
			for (int len = 1; len <= 35; len++) {
				for (int j = 0; j < LEN; j++) {
					to_calc.push_back(seq[j]);
				}
				double theta = kernel.NewThetaOptimizer(to_calc, Theta);
				double f = kernel.Fidelity(to_calc, theta);
				pair<double, double> cur = { fabs(theta - NeededAngle), f };
				if (cur < res) {
					res = cur;
					bestSeq = seq;
					bestLen = len;
					bestTheta = theta;
					bestF = f;
				}
			}
			return;
		}
		seq[pos] = -1;
		calc(pos + 1);

		seq[pos] = 0;
		calc(pos + 1);

		seq[pos] = 1;
		calc(pos + 1);
	};
	calc(0);
	kernel.WriteSequence(bestSeq);
	cout << bestLen << '\n';
	cout << fixed << setprecision(8) << "Theta = " << bestTheta << '\n' << "F = " << bestF << '\n';
}

void genetic(int CellsNumber) {
	int maxIter = 500;
	double crossover_probability = 0.8;
	double mutation_probability = 0.8;

	const double Len = (double)CellsNumber * T * NumberOfCycles;
	auto Amps = kernel.CreateAmpThresholds(CellsNumber);
	auto seq = kernel.CreateStartSCALLOP(CellsNumber, Amps[Amps.size() / 2]);
	CalculationDescriptor desc(w01, w12, wt, w, T, Len, tstep, Theta, NeededAngle, NumberOfCycles);
	GeneticAlgorithm<TYPE> algo(seq, crossover_probability, mutation_probability, maxIter, desc);
	algo.run();
	auto A1Sequence = algo.getSequence();
	auto A1F = algo.getLeak();
	auto A1Angle = algo.getAngle();
	auto nos = algo.getNumberOfCycles();
	cout << "N = " << CellsNumber << "; Theta = " << fixed << setprecision(10) << A1Angle << "; F = " << A1F << '\n';
	kernel.WriteSequence(A1Sequence);
}

int main() {
	setlocale(LC_ALL, "Rus");
	cout.precision(20);
	genetic(120);
	
	return 0;
}