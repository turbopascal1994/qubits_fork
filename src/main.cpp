#include <algorithm>
#include <cmath>
#include <chrono>
#include <complex>
#include <iostream>
#include <functional>
#include <fstream>
#include <limits.h>
#include <omp.h>
#include <map>
#include <vector>
#include "GeneticAlgorithm.h"

using namespace std;

void BruteForce(int CellsNumber = 120, int MaxCells = 120,
				double w01Coeff = 4, double NeededAngle = 0.032,
				string StringType = "bipolar") {
	int Type = (StringType == "bipolar" ? 3 : 2);
	double w01 = w01Coeff * 2 * PI * 1e9;
	const double w12 = w01 - 0.25 * 2 * PI * 1e9;
	const double wt = 25 * 2 * PI * 1e9;
	const double w = 4e-12;
	const double T = PI / wt * 2;
	const double Theta = 0.001;
	const double tstep = 5e-14;

	const int NumberOfCycles = (MaxCells + CellsNumber - 1) / CellsNumber;

	Kernel kernel(tstep, w01, w12, wt, w, T, Type);
	vector<int> seq(CellsNumber);
	pair<double, double> res = { INT_MAX, INT_MAX };
	vector<int> bestSeq;
	int bestLen;
	double bestTheta, bestF;
	function<void(int)> calc = [&](int pos) {
		if (pos == CellsNumber) {
			vector<int> to_calc;
			for (int len = 1; len <= NumberOfCycles; ++len) {
				for (int j = 0; j < CellsNumber; ++j) {
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
	kernel.WriteSequence(cout, bestSeq, '\t');
	cout << bestLen << '\t';
	cout << fixed << setprecision(8) << "Theta = " << bestTheta << '\t' << "F = " << bestF << '\n';
}

void Genetic(int CellsNumber=120, int MaxCells=120,
			 double w01Coeff=3,
			 double w12Coeff=0.25,
			 double wtCoeff = 25,
			 double NeededAngle=0.024,
			 double FidelityUpperBound=1e-4,
			 double CrossoverProbability=0.8,
			 double MutationProbability=0.8,
			 int MaxIter=500,
			 string StringType="bipolar") {
	int Type = (StringType == "bipolar" ? 3 : 2);
	double w01 = w01Coeff * 2 * PI * 1e9;
	double w12 = w01 - w12Coeff * 2 * PI * 1e9;
	double wt = wtCoeff * 2 * PI * 1e9;
	const double w = 4e-12;
	const double T = PI / wt * 2;
	const double Theta = 0.001;
	const double tstep = 5e-14;

	const int NumberOfCycles = (MaxCells + CellsNumber - 1) / CellsNumber;

	Kernel kernel(tstep, w01, w12, wt, w, T, Type);

	auto Amps = kernel.CreateAmpThresholds(CellsNumber);
	auto seq = kernel.CreateStartSCALLOP(CellsNumber, Amps[Amps.size() / 2]);
	CalculationDescriptor desc(w01, w12, wt, w, T, tstep, Theta, NeededAngle, NumberOfCycles, FidelityUpperBound, Type);
	GeneticAlgorithm algo(seq, CrossoverProbability, MutationProbability, MaxIter, desc);
	auto exec_time = algo.run();
	auto A1Sequence = algo.getSequence();
	auto A1F = algo.getLeak();
	auto A1Angle = algo.getAngle();
	auto nos = algo.getNumberOfCycles();

	string filename = "result_" + to_string(CellsNumber) + "_" + to_string(MaxCells) + "_" + to_string(w01Coeff) +
						"_" + to_string(w12Coeff) + "_" + to_string(wtCoeff) +
						"_" + to_string(NeededAngle) + "_" + to_string(FidelityUpperBound) + ".txt";
	ofstream fout;
	fout.open(filename, std::ios::app);

	fout << CellsNumber << '\t';
	kernel.WriteSequence(fout, A1Sequence, '\t');
	fout << nos << '\t';
	fout << fixed << setprecision(20) << A1Angle << '\t' << A1F << '\t' << exec_time << '\n';
	fout.close();
}

int main(int argc, char** argv) {
	int CellsNumber = atoi(argv[1]);
	int MaxCells = atoi(argv[2]);
	double w01Coeff = atof(argv[3]);
	double w12Coeff = atof(argv[4]);
	double wtCoeff = atof(argv[5]);
	double NeededAngle = atof(argv[6]);
	double FidelityUpperBound = atof(argv[7]);
	double CrossoverProbability = atof(argv[8]);
	double MutationProbability = atof(argv[9]);
	int MaxIter = atoi(argv[10]);
	string StringType = argv[11];
	Genetic(
		CellsNumber,
		MaxCells,
		w01Coeff,
		w12Coeff,
		wtCoeff,
		NeededAngle,
		FidelityUpperBound,
		CrossoverProbability,
		MutationProbability,
		MaxIter,
		StringType
	);
	
	return 0;
}