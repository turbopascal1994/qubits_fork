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
#include <numeric>

using namespace std;

void Genetic(int CellsNumber = 120, int MaxCells = 120,
	double w01Coeff = 3,
	double w12Coeff = 0.25,
	double wtCoeff = 25,
	double NeededAngle = 0.024,
	double AngleUpperBound = 1e-4,
	double CrossoverProbability = 0.8,
	double MutationProbability = 0.8,
	int MaxIter = 500,
	int Type = 3) {
	cout << "CellsNumber = " << CellsNumber << endl;
	cout << "MaxCells = " << MaxCells << endl;
	cout << "w01Coeff = " << w01Coeff << endl;
	cout << "w12Coeff = " << w12Coeff << endl;
	cout << "wtCoeff = " << wtCoeff << endl;
	cout << "NeededAngle = " << NeededAngle << endl;
	cout << "AngleUpperBound = " << AngleUpperBound << endl;
	cout << "CrossoverProbability = " << CrossoverProbability << endl;
	cout << "MutationProbability = " << MutationProbability << endl;
	cout << "MaxIter = " << MaxIter << endl;
	cout << "Type = " << Type << endl;
	double w01 = w01Coeff * 2 * PI * 1e9;
	double w12 = w01 - w12Coeff * 2 * PI * 1e9;
	double wt = wtCoeff * 2 * PI * 1e9;
	const double w = 4e-12;
	const double T = PI / wt * 2;
	const double Theta = 0.001;
	const double tstep = 5e-14;

	const int NumberOfCycles = (MaxCells + CellsNumber - 1) / CellsNumber;

	Kernel kernel(tstep, w01, w12, wt, w, T);

	auto Amps = kernel.CreateAmpThresholds(CellsNumber);
	auto seq = kernel.CreateStartSCALLOP(CellsNumber, Amps[Amps.size() / 2], Type);
	CalculationDescriptor desc(w01, w12, wt, w, T, tstep, Theta, NeededAngle, NumberOfCycles, AngleUpperBound, Type);
	GeneticAlgorithm algo(seq, CrossoverProbability, MutationProbability, MaxIter, desc);
	auto exec_time = algo.run();
	auto A1Sequence = algo.getSequence();
	auto A1F = algo.getLeak();
	auto A1Angle = algo.getAngle();
	auto nos = algo.getNumberOfCycles();

	string filename = "L=" + to_string(CellsNumber) + "_MaxL=" + to_string(MaxCells) + "_w01=" + to_string(w01Coeff) +
		"_w12=" + to_string(w12Coeff) + "_wt=" + to_string(wtCoeff) +
		"_Angle=" + to_string(NeededAngle) + "_" + to_string(AngleUpperBound) + ".txt";
	ofstream fout;
	fout.open(filename, std::ios::app);

	fout << CellsNumber << '\t';
	kernel.WriteSequence(fout, A1Sequence, '\t');
	fout << nos << '\t';
	fout << fixed << setprecision(20) << A1Angle << '\t';
	fout << A1F.leak << ' ' << 1 - A1F.fidelity << '\t';
	fout << '\t' << exec_time << '\n';
	fout.close();
}

map<string, double> preproc_args(int argc, char** argv) {
	map<string, double> mp = {
	  {"len", 120},
	  {"max_len", 120},
	  {"w01", 3},
	  {"w12", 0.25},
	  {"wt", 25},
	  {"angle", 0.024},
	  {"module", 0.0001},
	  {"mp", 0.8},
	  {"cp", 0.8},
	  {"iter", 500},
	  {"type", 3}
	};
	for (int i = 1; i < argc; i += 2) {
		string name = argv[i];
		name = name.substr(2, string::npos);
		double value;
		if (name == "type") {
			string val = argv[i + 1];
			value = val == "bipolar" ? 3 : 2;
		}
		else {
			value = atof(argv[i + 1]);
		}
		mp[name] = value;
	}
	return mp;
}

int main(int argc, char** argv) {
	omp_set_num_threads(4);
	auto mp = preproc_args(argc, argv);
	Genetic(
		int(mp["len"]),
		int(mp["max_len"]),
		mp["w01"],
		mp["w12"],
		mp["wt"],
		mp["angle"],
		mp["module"],
		mp["cp"],
		mp["mp"],
		mp["iter"],
		int(mp["type"])
	);
	return 0;
}