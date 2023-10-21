#include "GeneticAlgorithm.h"
#include "ArgsPreprocessor.h"

using namespace std;

void Genetic(
	int RegularLen = 30,
	int CellsNumber = 120, 
	int MaxCells = 120,
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

	const int NumberOfCycles = 1;

	ConstantsDescriptor config(w01, w12, wt, w, T, tstep, Theta, NeededAngle, NumberOfCycles, AngleUpperBound, Type, RegularLen);

	vector<vector<int>> seqs(2 * CellsNumber);
	uniform_int_distribution<> dist(-1, 1);
	random_device rd;
	mt19937 gen(rd());
	for (auto& seq : seqs) {
		seq.resize(CellsNumber);
		for (auto& j : seq) {
			j = dist(gen);
		}
	}
	GeneticHyperParameters hyperParams(CrossoverProbability, MutationProbability, MaxIter);
	GeneticAlgorithm algo(seqs, config, hyperParams);
	auto exec_time = algo.run();

	string filename = "RL=" + to_string(RegularLen) + "_L=" + to_string(CellsNumber) + "_w01=" + to_string(w01Coeff) +
		"_w12=" + to_string(w12Coeff) + "_wt=" + to_string(wtCoeff) +
		"_Angle=" + to_string(NeededAngle) + ".txt";
	ofstream fout;
	fout.open(filename, std::ios::app);

	fout << CellsNumber << '\t';
	for (auto& i : algo.getSequence()) fout << i;
	fout << '\t';
	fout << algo.getNumberOfCycles() << '\t';
	fout << algo.getBestIteration() << '\t';
	fout << 1 - algo.getFidelity()  << ' ' << algo.getLeak() << '\t';
	fout << '\t' << exec_time << '\n';
	fout.close();
}

int main(int argc, char** argv) {
	omp_set_num_threads(4);
	auto mp = ArgsPreprocessor::run(argc, argv);
	Genetic(
		int(mp["regular_len"]),
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