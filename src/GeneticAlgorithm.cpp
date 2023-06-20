#include "GeneticAlgorithm.h"

GeneticAlgorithm::GeneticAlgorithm(
	const std::vector<std::vector<int>>& _sequences, 
	ConstantsDescriptor _config, 
	GeneticHyperParameters _hyperParams
) : config(_config), kernel(_config) {
	hyperParams = _hyperParams;
	populationSize = _sequences.size();
	for (size_t i = 0; i < _sequences.size(); ++i) {
		population.push_back(CreateIndividual(_sequences[i]));
	}
}

double GeneticAlgorithm::getLeak() {
	return population[0].leak;
}

double GeneticAlgorithm::getFidelity() {
	return population[0].fidelity;
}

vector<int> GeneticAlgorithm::getSequence() {
	return population[0].sequence;
}

int GeneticAlgorithm::getNumberOfCycles() {
	return population[0].numberOfCycles;
}

BaseIndividual GeneticAlgorithm::CreateIndividual(const std::vector<int>& sequence) {
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

Kernel::FidelityResult GeneticAlgorithm::_compute_fidelity(std::vector<int>& sequence, double neededAngle) {
	return kernel.Fidelity(sequence, neededAngle);
}

void GeneticAlgorithm::CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs) {
	assert(ls.size() == rs.size());
	int index = GenerateInt(2, (int)rs.size() - 1);
	for (size_t i = index; i < ls.size(); ++i) {
		swap(ls[i], rs[i]);
	}
}

void GeneticAlgorithm::MutationImpl(std::vector<int>& sequence) {
	double p = 1.0 / sequence.size();
	vector<int> tmp;
	for (int i = 0; i < sequence.size(); ++i) {
		if (GenerateProbability() < p) {
			if (config.type == 2) tmp = { 0, 1 };
			else tmp = { -1, 0, 1 };
			tmp.erase(find(tmp.begin(), tmp.end(), sequence[i]));
			if (tmp.size() == 1) {
				sequence[i] = tmp[0];
			}
			else {
				sequence[i] = tmp[GenerateProbability() * 2 < 1.0];
			}
		}
	}
}