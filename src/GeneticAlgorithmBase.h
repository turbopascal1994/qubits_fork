#pragma once
#include<vector>
#include<array>
#include<algorithm>
#include<limits>
#include<random>
#include<omp.h>
#include "Individuals.h"

struct GeneticHyperParameters {
	double crossoverProbability;
	double mutationProbability;
	unsigned int maxIter;
	GeneticHyperParameters(
		double _crossoverProbability = 0.8,
		double _mutationProbability = 0.8,
		unsigned int _maxIter = 500
	) : crossoverProbability(_crossoverProbability),
		mutationProbability(_mutationProbability),
		maxIter(_maxIter) {}
};

template<class INDIVID>
class GeneticAlgorithmBase {
protected:
	std::vector<INDIVID> population;
	unsigned int populationSize;
	virtual INDIVID CreateIndividual(const std::vector<int>& sequence) = 0;
	GeneticHyperParameters hyperParams;


	// START OF RANDOM BLOCK
	std::mt19937 randomGenerator;
	int GenerateInt(int l, int r) {
		std::uniform_int_distribution<> dist(l, r);
		return dist(randomGenerator);
	}
	double GenerateProbability() {
		std::uniform_int_distribution<> dist(0, std::numeric_limits<int>::max());
		return 1.0 * dist(randomGenerator) / std::numeric_limits<int>::max();
	}
	// END OF RANDOM BLOCK

	std::array<std::vector<int>, 2> Tournament() {
		int index1, index2, index3;
		do {
			index1 = GenerateInt(0, populationSize - 1);
			index2 = GenerateInt(0, populationSize - 1);
			index3 = GenerateInt(0, populationSize - 1);
		} while (index1 == index2 || index1 == index3 || index2 == index3);
		vector<INDIVID> indices = { 
			population[index1], 
			population[index2], 
			population[index3] 
		};
		std::sort(indices.begin(), indices.end());
		return { indices[0].sequence, indices[1].sequence };
	}

	// START OF CROSSOVER BLOCK
	void Crossover(std::vector<int>& ls, std::vector<int>& rs) {
		CrossoverImpl(ls, rs);
	}
	virtual void CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs) = 0;
	// END OF CROSSOVER BLOCK

	// START OF MUTATION BLOCK
	void Mutation(std::vector<int>& sequence) {
		MutationImpl(sequence);
	}
	virtual void MutationImpl(std::vector<int>& sequence) = 0;
	// END OF MUTATION BLOCK

public:
	GeneticAlgorithmBase() : randomGenerator(chrono::high_resolution_clock::now().time_since_epoch().count()) {}
	double run() {
		double start = omp_get_wtime();
		int offSpringSize = populationSize / 2;
		std::vector<INDIVID> newPopulation(populationSize);

		for (int iteration = 0; iteration < hyperParams.maxIter; ++iteration) {
#pragma omp parallel for shared(newPopulation)
			for (int i = 0; i < offSpringSize; ++i) {
				auto parents = Tournament();
				if (GenerateProbability() < hyperParams.crossoverProbability) {
					Crossover(parents[0], parents[1]);
				}
				for (int j = 0; j < 2; ++j) {
					if (GenerateProbability() < hyperParams.mutationProbability) {
						Mutation(parents[j]);
					}
					newPopulation[i * 2 + j] = CreateIndividual(parents[j]);
				}
			}
			std::sort(population.begin(), population.end());
			std::sort(newPopulation.begin(), newPopulation.end());
			newPopulation[populationSize - 2] = population[0];
			newPopulation[populationSize - 1] = population[1];
			swap(population, newPopulation);
		}
		sort(population.begin(), population.end());
		return omp_get_wtime() - start;
	}
};