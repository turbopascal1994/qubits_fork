#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <cassert>
#include <omp.h>
#include "kernel.h"
#include <limits.h>
#include "GeneticAlgorithmBase.h"
#include "ConstantsDescriptor.h"

class GeneticAlgorithm : public GeneticAlgorithmBase<BaseIndividual> {
public:
	GeneticAlgorithm(
		const std::vector<std::vector<int>>& _sequences,
		ConstantsDescriptor _config,
		GeneticHyperParameters _hyperParams
	);
	double getLeak();
	double getFidelity();
	vector<int> getSequence();
	int getNumberOfCycles();
private:
	ConstantsDescriptor config;
	Kernel kernel;

	BaseIndividual CreateIndividual(const std::vector<int>& sequence);
	Kernel::FidelityResult _compute_fidelity(std::vector<int>& sequence, double neededAngle);
	void CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs);
	void MutationImpl(std::vector<int>& sequence);
	bool CheckStopCondition();
	vector<int> buildSequence(const vector<int>& initial, int prefix_count);
};
