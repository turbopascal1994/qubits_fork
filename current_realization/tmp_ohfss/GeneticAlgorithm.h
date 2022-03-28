#pragma once
#include <vector>
#include <random>
#include <chrono>
#include "kernel.h"

using namespace std;

struct Individual {
	vector<int> sequence;
	double optimizedTheta;
	double F;
};

struct CalculationDescriptor {
	double w01, w12, w, T, Len, tstep, Theta;
	int numberOfCycles;
	double neededAngle;
};

class GeneticAlgorithm
{
	vector<int> inputSequence;
	double crossover_probability;
	double mutation_probability;
	int maxIter;
	int population_size;
	vector<Individual> population;
	CalculationDescriptor config;
	mt19937 gen;
	Kernel kernel;
	double bestAbs, bestLeak, bestAngle;
	vector<int> bestSequence;

	void createPopulation(vector<int> sequence);
	Individual createIndividual(vector<int> sequence);
	double _optimizedTheta(vector<int>& Sequence);
	double _fidelity(vector<int>& sequence, double theta);
	vector<Individual> tournament();
	double generateProbability();
	int generateInt(int l, int r);
	void crossover(Individual& a, Individual& b);
	void mutation(Individual& a, double p);
public:
	GeneticAlgorithm(std::vector<int> _inputSequence, double _crossover_probability,
					 double _mutation_probability, int _maxIter, CalculationDescriptor _config);
	void run();
	void reset();
	double getLeak();
	double getAngle();
	vector<int> getSequence();
};
