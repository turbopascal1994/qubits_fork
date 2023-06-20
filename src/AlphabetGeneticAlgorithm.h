#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <cassert>
#include <omp.h>
#include "kernel.h"
#include <limits.h>
#include "GeneticAlgorithmBase.h"
#include "AlphabetConstantsDescriptor.h"

class AlphabetGeneticAlgorithm : public GeneticAlgorithmBase<BaseIndividual> {
public:
	AlphabetGeneticAlgorithm(
		const std::vector<std::vector<int>>& _sequences,
		AlphabetConstantsDescriptor _config,
		GeneticHyperParameters _hyperParams
	);
	double getLeak();
	double getFidelity();
	vector<int> getSequence();
	int getNumberOfCycles();
private:
	AlphabetConstantsDescriptor config;
	Kernel kernel;

	BaseIndividual CreateIndividual(const std::vector<int>& sequence);
	Kernel::FidelityResult _compute_fidelity(std::vector<int>& sequence, double neededAngle);
	void CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs);
	void MutationImpl(std::vector<int>& sequence);
};
