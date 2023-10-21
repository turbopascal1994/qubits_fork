#include "Alphabet.h"
#include "AlphabetGeneticAlgorithm.h"
#include "ArgsPreprocessor.h"

using namespace std;

void Genetic(
	int CellsNumber,
	double w01Coeff,
	double w12Coeff,
	double wtCoeff,
	double NeededAngle,
	double AngleUpperBound,
	double CrossoverProbability,
	double MutationProbability,
	int MaxIter,
	int Type
) {
	// cout << "CellsNumber = " << CellsNumber << endl;
	// cout << "w01Coeff = " << w01Coeff << endl;
	// cout << "w12Coeff = " << w12Coeff << endl;
	// cout << "wtCoeff = " << wtCoeff << endl;
	// cout << "NeededAngle = " << NeededAngle << endl;
	// cout << "AngleUpperBound = " << AngleUpperBound << endl;
	// cout << "CrossoverProbability = " << CrossoverProbability << endl;
	// cout << "MutationProbability = " << MutationProbability << endl;
	// cout << "MaxIter = " << MaxIter << endl;
	// cout << "Type = " << Type << endl;
	double w01 = w01Coeff * 2 * PI * 1e9;
	double w12 = w01 - w12Coeff * 2 * PI * 1e9;
	double wt = wtCoeff * 2 * PI * 1e9;
	const double T = PI / wt * 2;
	const double Theta = 0.001;
	const double tstep = 5e-14;
	const double w = tstep;

	int populationSize = 4 * CellsNumber;
	int alphabetSize = 100;
	Alphabet alphabet(alphabetSize, Type);
	random_device rd;
	mt19937 gen(rd());

	string filename = "L=" + to_string(CellsNumber) + "_w01=" + to_string(w01Coeff) +
		"_w12=" + to_string(w12Coeff) + "_wt=" + to_string(wtCoeff) + ".txt";
	ofstream fout;
	fout.open(filename, std::ios::app);

	fout << "ALPHABET START\n";
	fout << alphabet.wordbook.size() << '\n';
	for (auto & i : alphabet.wordbook) {
		for(int j = 0;j < 3;j++){
			fout << i[j] << ' ';
		}
		fout << '\n';
	}
	fout << "ALPHABET END\n";

	for (int NumberOfCycles = 1; NumberOfCycles <= 10; NumberOfCycles++) {

		vector<vector<int>> seqs(populationSize);
		uniform_int_distribution<> dist(0, alphabetSize - 1);
		for (auto& seq : seqs) {
			seq.resize(CellsNumber);
			for (auto& j : seq) {
				j = dist(gen);
			}
		}
		AlphabetConstantsDescriptor config(alphabet, w01, w12, wt, w, T, tstep, Theta, NeededAngle, NumberOfCycles, AngleUpperBound, Type);
		GeneticHyperParameters hyperParams(CrossoverProbability, MutationProbability, MaxIter);
		AlphabetGeneticAlgorithm algo(seqs, config, hyperParams);
		auto exec_time = algo.run();

		auto sequence = algo.getSequence();
		fout << sequence.size() << '\t';
		for (auto& i : sequence) fout << i << ' ';
		fout << '\t';
		auto decoded_sequence = alphabet.decode(sequence);
		fout << decoded_sequence.size() << '\t';
		for (auto& i : decoded_sequence) fout << i;
		fout << '\t';
		fout << algo.getNumberOfCycles() << '\t';
		fout << algo.getBestIteration() << '\t';
		fout << 1 - algo.getFidelity() << ' ' << algo.getLeak() << '\t';
		fout << '\t' << exec_time << '\n';
	}
	fout.close();
}

int main(int argc, char** argv) {
	omp_set_num_threads(4);
	auto mp = ArgsPreprocessor::run(argc, argv);
	Genetic(
		int(mp["len"]),
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