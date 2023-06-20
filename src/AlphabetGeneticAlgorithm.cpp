#include "AlphabetGeneticAlgorithm.h"
#include <numeric>

AlphabetGeneticAlgorithm::AlphabetGeneticAlgorithm(
	const std::vector<std::vector<int>>& _sequences,
	AlphabetConstantsDescriptor _config,
	GeneticHyperParameters _hyperParams
) : config(_config), kernel(_config) {
	hyperParams = _hyperParams;
	populationSize = _sequences.size();
	for (size_t i = 0; i < _sequences.size(); ++i) {
		population.push_back(CreateIndividual(_sequences[i]));
	}
}

double AlphabetGeneticAlgorithm::getLeak() {
	return population[0].leak;
}

double AlphabetGeneticAlgorithm::getFidelity() {
	return population[0].fidelity;
}

vector<int> AlphabetGeneticAlgorithm::getSequence() {
	return population[0].sequence;
}

int AlphabetGeneticAlgorithm::getNumberOfCycles() {
	return population[0].numberOfCycles;
}

BaseIndividual AlphabetGeneticAlgorithm::CreateIndividual(const std::vector<int>& sequence) {
	vector<int> cur_seq;
	cur_seq.reserve(sequence.size() * config.numberOfCycles);

	double fidelity = -1;
	double leak = -1;
	int NumberOfCycles = -1;

	for (int len = 1; len <= config.numberOfCycles; ++len) {
		for (size_t j = 0; j < sequence.size(); ++j) {
			cur_seq.push_back(sequence[j]);
		}

		Kernel::FidelityResult res = _compute_fidelity(cur_seq, config.neededAngle);
		double cur_fidelity = res.fidelity;
		if (cur_fidelity > fidelity || fidelity == -1) {
			fidelity = cur_fidelity;
			leak = res.leak;
			NumberOfCycles = len;
		}
	}
	return BaseIndividual(sequence, fidelity, leak, NumberOfCycles);
}

Kernel::FidelityResult AlphabetGeneticAlgorithm::_compute_fidelity(std::vector<int>& sequence, double neededAngle) {
	return kernel.Fidelity(config.alp.decode(sequence), neededAngle);
}

void AlphabetGeneticAlgorithm::CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs) {
	assert(ls.size() == rs.size());
	int index = GenerateInt(2, (int)rs.size() - 1);
	for (size_t i = index; i < ls.size(); ++i) {
		swap(ls[i], rs[i]);
	}
}

void AlphabetGeneticAlgorithm::MutationImpl(std::vector<int>& sequence) {
	double p = 1.0 / sequence.size();
	int alphabetSize = config.alp.wordbook.size();
	vector<int> tmp(alphabetSize);
	std::iota(tmp.begin(), tmp.end(), 0);
	for (int i = 0; i < sequence.size(); ++i) {
		if (GenerateProbability() < p) {
			do {
				shuffle(tmp.begin(), tmp.end(), randomGenerator);
			} while (tmp[0] == sequence[i]);
			sequence[i] = tmp[0];
		}
	}
}
