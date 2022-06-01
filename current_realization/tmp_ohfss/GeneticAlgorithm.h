#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <omp.h>
#include "kernel.h"

using namespace std;

struct Individual {
	vector<int> sequence;
	double optimizedTheta;
	double F;
	int numberOfCycles;
	Individual(vector<int> _sequence = {0}, double _theta = 0, double _f = 0, int _cycles = 0) :
		sequence(_sequence), optimizedTheta(_theta), F(_f), numberOfCycles(_cycles) {}
};

struct CalculationDescriptor {
	double w01, w12, wt, w, T, Len, tstep, Theta;
	int numberOfCycles;
	double neededAngle;
	CalculationDescriptor(
		double w01, double w12, double wt, double w,
		double T, double Len, double tstep, double Theta,
		double neededAngle, int numberOfCycles) :
		w01(w01), w12(w12), wt(wt), w(w),
		T(T), Len(Len), tstep(tstep), Theta(Theta),
		neededAngle(neededAngle), numberOfCycles(numberOfCycles) {}
};

template<int TYPE = 3>
class GeneticAlgorithm
{
	double crossover_probability;
	double mutation_probability;
	int maxIter;
	int population_size;
	vector<Individual> population;
	CalculationDescriptor config;
	mt19937 gen;
	Kernel<TYPE> kernel;
	double bestLeak, bestAngle;
	int bestNumberOfCycles;
	vector<int> bestSequence;

	void createPopulation(vector<int> sequence);
	Individual createIndividual(vector<int> sequence);
	double _optimizedTheta(vector<int>& Sequence);
	double _fidelity(vector<int>& sequence, double theta);
	vector<vector<int>> tournament();
	double generateProbability();
	int generateInt(int l, int r);
	void crossover(vector<int>& a, vector<int>& b);
	void mutation(vector<int>& a, double p);
	int compare(Individual& a, Individual& b);
	void sort(vector<Individual>& a);
public:
	GeneticAlgorithm(std::vector<int> _inputSequence, double _crossover_probability,
					 double _mutation_probability, int _maxIter, CalculationDescriptor _config);
	GeneticAlgorithm(std::vector<vector<int>> _inputPopulation, double _crossover_probability,
					 double _mutation_probability, int _maxIter, CalculationDescriptor _config);
	void run();
	void reset();
	double getLeak();
	double getAngle();
	vector<int> getSequence();
	int getNumberOfCycles();
};

template<int TYPE>
GeneticAlgorithm<TYPE>::GeneticAlgorithm(
	std::vector<int> _inputSequence,
	double _crossover_probability,
	double _mutation_probability,
	int _maxIter,
	CalculationDescriptor _config
) :
	crossover_probability(_crossover_probability),
	mutation_probability(_mutation_probability), maxIter(_maxIter), config(_config),
	kernel(config.tstep, config.w01, config.w12, config.wt, config.w, config.T),
	gen(mt19937(chrono::high_resolution_clock::now().time_since_epoch().count())) {
	reset();
	createPopulation(_inputSequence);
}

template<int TYPE>
GeneticAlgorithm<TYPE>::GeneticAlgorithm(
	std::vector<vector<int>> _inputPopulation,
	double _crossover_probability,
	double _mutation_probability,
	int _maxIter,
	CalculationDescriptor _config
) :
	crossover_probability(_crossover_probability),
	mutation_probability(_mutation_probability), maxIter(_maxIter), config(_config),
	kernel(config.tstep, config.w01, config.w12, config.wt, config.w, config.T),
	gen(mt19937(chrono::high_resolution_clock::now().time_since_epoch().count())) {
	reset();
	population_size = _inputPopulation.size();
	for (int i = 0; i < _inputPopulation.size(); i++) {
		population.push_back(createIndividual(_inputPopulation[i]));
	}
}

template<int TYPE>
void GeneticAlgorithm<TYPE>::createPopulation(vector<int> sequence) {
	population.clear();
	population.push_back(createIndividual(sequence));
	for (int i = 0; i < sequence.size() * 2; i++) {
		population.push_back(createIndividual(kernel.ChangingOneElement(sequence, i / 2, i % 2)));
	}
	population_size = population.size();
}

template<int TYPE>
Individual GeneticAlgorithm<TYPE>::createIndividual(vector<int> sequence) {
	vector<int> cur_seq;
	cur_seq.reserve(sequence.size() * config.numberOfCycles);

	pair<double, double> Q = {INT_MAX, INT_MAX};
	int numberOfCycles = -1;
	double theta = INT_MAX;

	for (int len = 1; len <= config.numberOfCycles; len++) {
		for (size_t j = 0; j < sequence.size(); j++) cur_seq.push_back(sequence[j]);
		double optimizedTheta = _optimizedTheta(cur_seq);
		double F = _fidelity(cur_seq, optimizedTheta);
		pair<double, double> cur_Q = { fabs(optimizedTheta - config.neededAngle), F };
		if (cur_Q < Q) {
			Q = cur_Q;
			theta = optimizedTheta;
			numberOfCycles = len;
		}
	}
	return Individual(sequence, theta, Q.second, numberOfCycles);
}

template<int TYPE>
int GeneticAlgorithm<TYPE>::compare(Individual& a, Individual& b) {
	double left_abs = fabs(a.optimizedTheta - config.neededAngle);
	double right_abs = fabs(b.optimizedTheta - config.neededAngle);
	if (fabs(left_abs - right_abs) < 1e-9) {
		if (a.F < b.F) return -1;
		return 1;
	}
	if (left_abs < right_abs) return -1;
	return 1;
}

template<int TYPE>
void GeneticAlgorithm<TYPE>::sort(vector<Individual> & a) {
	std::sort(a.begin(), a.end(), [&](Individual& l, Individual& r) {
		return compare(l, r) == -1;
	});
}

template<int TYPE>
vector<vector<int>> GeneticAlgorithm<TYPE>::tournament() {
	int index1, index2, index3;
	do {
		index1 = generateInt(0, population_size - 1);
		index2 = generateInt(0, population_size - 1);
		index3 = generateInt(0, population_size - 1);
	} while (index1 == index2 || index1 == index3 || index2 == index3);
	vector<vector<int>> ans;
	if (compare(population[index1], population[index2]) == -1) ans.push_back(population[index1].sequence);
	else ans.push_back(population[index2].sequence);

	if (compare(population[index2], population[index3]) == -1) ans.push_back(population[index2].sequence);
	else ans.push_back(population[index3].sequence);
	return ans;
}

template<int TYPE>
double GeneticAlgorithm<TYPE>::generateProbability(){
	double g = gen();
	g = fabs(g);
	return g / numeric_limits<int>::max();
}

template<int TYPE>
int GeneticAlgorithm<TYPE>::generateInt(int l, int r) {
	int q = gen();
	q = abs(q);
	return l + q % (r - l + 1);
}

template<int TYPE>
void GeneticAlgorithm<TYPE>::crossover(vector<int>& a, vector<int>& b) {
	int index = generateInt(2, (int)b.size() - 3);
	for (int i = index; i < a.size(); i++) {
		swap(a[i], b[i]);
	}
}

template<int TYPE>
void GeneticAlgorithm<TYPE>::mutation(vector<int>& a, double p) {
	for (int i = 0; i < a.size(); i++) {
		if (generateProbability() < p) {
			vector<int> tmp;
			if (TYPE == 2) tmp = { 0, 1 };
			else tmp = { -1, 0, 1 };
			tmp.erase(find(tmp.begin(), tmp.end(), a[i]));
			if (tmp.size() == 1) a[i] = tmp[0];
			else a[i] = tmp[generateProbability() < 1.0 / 2];
		}
	}
}

template<int TYPE>
void GeneticAlgorithm<TYPE>::run() {
	double start = omp_get_wtime();
	vector<vector<int>> offSpring;
	offSpring.reserve(population_size);

	for (int iteration = 0; iteration < maxIter; iteration++) {
		// cout << iteration << '\n';
		offSpring.clear();
		for (int i = 0; i < population_size; i += 2) {
			auto parents = tournament();
			if (generateProbability() < crossover_probability) {
				crossover(parents[0], parents[1]);
			}
			for (int j = 0; j < 2; j++) {
				if (generateProbability() < mutation_probability) {
					mutation(parents[j], 1.0 / population[0].sequence.size());
				}
				offSpring.push_back(parents[j]);
			}
		}
		vector<Individual> new_population(offSpring.size());
#pragma omp parallel for
		for (int i = 0; i < new_population.size(); i++) {
			new_population[i] = createIndividual(offSpring[i]);
		}
		sort(population);
		sort(new_population);
		new_population[population_size - 2] = population[0];
		new_population[population_size - 1] = population[1];
		population = move(new_population);
	}
	sort(population);
	bestLeak = population[0].F;
	bestAngle = population[0].optimizedTheta;
	bestSequence = population[0].sequence;
	bestNumberOfCycles = population[0].numberOfCycles;
	cout << "Время работы генетического алгоритма = " << omp_get_wtime() - start << '\n';
}

template<int TYPE>
double GeneticAlgorithm<TYPE>::_optimizedTheta(vector<int>& sequence)
{
	double res = kernel.NewThetaOptimizer(sequence, config.Theta);
	return res;
}

template<int TYPE>
double GeneticAlgorithm<TYPE>::_fidelity(vector<int>& sequence, double theta)
{
	double res = kernel.Fidelity(sequence, theta);
	return res;
}

template<int TYPE>
void GeneticAlgorithm<TYPE>::reset() {
	bestLeak = numeric_limits<double>::max();
	bestAngle = -1;
	bestSequence.clear();
}

template<int TYPE>
double GeneticAlgorithm<TYPE>::getLeak()
{
	return bestLeak;
}

template<int TYPE>
double GeneticAlgorithm<TYPE>::getAngle()
{
	return bestAngle;
}

template<int TYPE>
vector<int> GeneticAlgorithm<TYPE>::getSequence()
{
	return bestSequence;
}

template<int TYPE>
int GeneticAlgorithm<TYPE>::getNumberOfCycles()
{
	return bestNumberOfCycles;
}
