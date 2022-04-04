#include "GeneticAlgorithm.h"
#include <cassert>
#include<omp.h>

GeneticAlgorithm::GeneticAlgorithm(
	std::vector<int> _inputSequence,
	double _crossover_probability,
	double _mutation_probability,
	int _maxIter,
	CalculationDescriptor _config
):
		inputSequence(_inputSequence), crossover_probability(_crossover_probability),
		mutation_probability(_mutation_probability), maxIter(_maxIter), config(_config),
		kernel(config.tstep, config.w01, config.w12, config.wt, config.w, config.T){
	gen = mt19937(chrono::high_resolution_clock::now().time_since_epoch().count());
	population_size = 2 * inputSequence.size();
	reset();
	createPopulation(inputSequence);
}

void GeneticAlgorithm::createPopulation(vector<int> sequence){
	population.clear();
	for (int i = 0; i < population_size; i++) {
		population.push_back(createIndividual(kernel.ChangingOneElement(sequence, i / 2, i % 2)));
	}
}

Individual GeneticAlgorithm::createIndividual(vector<int> sequence){
	double optimizedTheta = _optimizedTheta(sequence);
	double F = _fidelity(sequence, optimizedTheta);
	return Individual({ sequence, optimizedTheta, F });
}

vector<Individual> GeneticAlgorithm::tournament(){
	vector<Individual> ans(population_size);
#pragma omp parallel for
	for (int i = 0; i < population_size; i++) {
		int index1, index2, index3;
		do {
			index1 = gen() % population_size;
			index2 = gen() % population_size;
			index3 = gen() % population_size;
		} while (index1 == index2 || index1 == index3 || index2 == index3);
		vector<pair<double, int>> to_sort;
		to_sort.push_back({ fabs(population[index1].optimizedTheta - config.neededAngle), index1 });
		to_sort.push_back({ fabs(population[index2].optimizedTheta - config.neededAngle), index2 });
		to_sort.push_back({ fabs(population[index3].optimizedTheta - config.neededAngle), index3 });
		sort(to_sort.begin(), to_sort.end());
		ans[i] = population[to_sort[0].second];
	}
	return ans;
}

double GeneticAlgorithm::generateProbability()
{
	double g = gen();
	return g / numeric_limits<int>::max();
}

int GeneticAlgorithm::generateInt(int l, int r){
	return l + gen() % (r - l + 1);
}

void GeneticAlgorithm::crossover(Individual& a, Individual& b){
	int index = generateInt(2, b.sequence.size() - 3);
#pragma omp parallel for
	for (int i = index; i < a.sequence.size(); i++) {
		swap(a.sequence[i], b.sequence[i]);
	}
	a = createIndividual(a.sequence);
	b = createIndividual(b.sequence);
}

void GeneticAlgorithm::mutation(Individual& a, double p){
	bool changed = false;
	for (int i = 0; i < a.sequence.size(); i++) {
		if (generateProbability() < p) {
			changed = true;
			vector<int> tmp = { -1, 0, 1 };
			tmp.erase(find(tmp.begin(), tmp.end(), a.sequence[i]));
			a.sequence[i] = tmp[generateProbability() < 1.0 / 2];
		}
	}
	if (changed) a = createIndividual(a.sequence);
}

void GeneticAlgorithm::run(){
	double start = omp_get_wtime();
	for (int iteration = 0; iteration < maxIter; iteration++) {
		auto offspring = tournament();
#pragma omp parallel for
		for (int i = 0; i < population_size; i += 2) {
			if (generateProbability() < crossover_probability) {
				crossover(offspring[i], offspring[i + 1]);
			}
		}
#pragma omp parallel for
		for (int i = 0; i < population_size; i++) {
			if (generateProbability() < mutation_probability) {
				mutation(offspring[i], 1.0 / inputSequence.size());
			}
		}
		for (int i = 0; i < population_size; i++) {
			if (bestAbs > fabs(offspring[i].optimizedTheta - config.neededAngle)) {
				bestAbs = fabs(offspring[i].optimizedTheta - config.neededAngle);
				bestLeak = offspring[i].F;
				bestAngle = offspring[i].optimizedTheta;
				bestSequence = offspring[i].sequence;
			}
		}
		population = offspring;
	}
	cout << "Run time = " << omp_get_wtime() - start << '\n';
}

double GeneticAlgorithm::_optimizedTheta(vector<int>& sequence)
{
	return kernel.NewThetaOptimizer(sequence, config.Theta);
}

double GeneticAlgorithm::_fidelity(vector<int>& sequence, double theta)
{
	return kernel.Fidelity(sequence, theta);
}

void GeneticAlgorithm::reset() {
	bestLeak = numeric_limits<double>::max();
	bestAbs = numeric_limits<double>::max();
	bestAngle = -1;
	bestSequence.clear();
}

double GeneticAlgorithm::getLeak()
{
	return bestLeak;
}

double GeneticAlgorithm::getAngle()
{
	return bestAngle;
}

vector<int> GeneticAlgorithm::getSequence()
{
	return bestSequence;
}
