#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <omp.h>
#include "kernel.h"
#include <limits.h>

using namespace std;

struct Individual {
	vector<int> Sequence;
	double OptimizedTheta;
	double Fidelity;
	int NumberOfCycles;
	Individual() {}
	Individual(const vector<int>& _sequence, double _theta, double _f, int _cycles) :
		Sequence(_sequence), OptimizedTheta(_theta), Fidelity(_f), NumberOfCycles(_cycles) {}
};

struct CalculationDescriptor {
	double w01, w12, wt, w, T, Tstep, Theta;
	int NumberOfCycles;
	double NeededAngle;
	double AngleUpperBound;
	int Type;
	CalculationDescriptor(
		double w01, double w12, double wt, double w,
		double T, double Tstep, double Theta,
		double NeededAngle, int NumberOfCycles, double AngleUpperBound, int Type) :
		w01(w01), w12(w12), wt(wt), w(w),
		T(T), Tstep(Tstep), Theta(Theta),
		NeededAngle(NeededAngle), NumberOfCycles(NumberOfCycles),
		AngleUpperBound(AngleUpperBound), Type(Type) {}
};

class GeneticAlgorithm
{
	double CrossoverProbability;
	double MutationProbability;
	int MaxIter;
	int PopulationSize;
	vector<Individual> Population;
	CalculationDescriptor Config;
	mt19937 RandomGenerator;
	Kernel Kernel;
	tuple<double, double, double, double> BestLeak;
	double BestAngle;
	int BestNumberOfCycles;
	vector<int> BestSequence;

	void _createPopulation(const vector<int>& sequence);
	Individual _createIndividual(const vector<int>& sequence);
	double _optimizedTheta(const vector<int>& Sequence);
	tuple<double, double, double, double> _fidelity(const vector<int>& sequence, double theta);
	vector<vector<int>> _tournament();
	double _generateProbability();
	int _generateInt(int l, int r);
	void _crossover(vector<int>& a, vector<int>& b);
	void _mutation(vector<int>& a, double p);
	bool _compare(const Individual& a, const Individual& b);
public:
	GeneticAlgorithm(const std::vector<int>& _inputSequence, double _crossover_probability,
					 double _mutation_probability, int _maxIter, CalculationDescriptor _config);
	GeneticAlgorithm(const std::vector<vector<int>>& _inputPopulation, double _crossover_probability,
					 double _mutation_probability, int _maxIter, CalculationDescriptor _config);
	double run();
	void reset();
	tuple<double, double, double, double> getLeak();
	double getAngle();
	vector<int> getSequence();
	int getNumberOfCycles();
};
