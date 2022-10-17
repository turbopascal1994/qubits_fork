#include "GeneticAlgorithm.h"

GeneticAlgorithm::GeneticAlgorithm(const std::vector<int>& _inputSequence,
								   double _crossover_probability,
								   double _mutation_probability,
								   int _maxIter,
								   CalculationDescriptor _config) :
	CrossoverProbability(_crossover_probability),
	MutationProbability(_mutation_probability), MaxIter(_maxIter), Config(_config),
	Kernel(Config.Tstep, Config.w01, Config.w12, Config.wt, Config.w, Config.T, Config.Type),
	RandomGenerator(mt19937(chrono::high_resolution_clock::now().time_since_epoch().count())) {
		reset();
		_createPopulation(_inputSequence);
}

GeneticAlgorithm::GeneticAlgorithm(const std::vector<vector<int>>& _inputPopulation,
								   double _crossover_probability,
								   double _mutation_probability,
								   int _maxIter,
								   CalculationDescriptor _config) :
	CrossoverProbability(_crossover_probability),
	MutationProbability(_mutation_probability), MaxIter(_maxIter), Config(_config),
	Kernel(Config.Tstep, Config.w01, Config.w12, Config.wt, Config.w, Config.T, Config.Type),
	RandomGenerator(mt19937(chrono::high_resolution_clock::now().time_since_epoch().count())) {
		reset();
		PopulationSize = _inputPopulation.size();
		for (int i = 0; i < _inputPopulation.size(); ++i) {
			Population.push_back(_createIndividual(_inputPopulation[i]));
		}
	}

double GeneticAlgorithm::run() {
	double start = omp_get_wtime();
	int offSpringSize = PopulationSize / 2;
	vector<Individual> NewPopulation(PopulationSize);

	for (int iteration = 0; iteration < MaxIter; ++iteration) {
#pragma omp parallel for shared(NewPopulation)
		for (int i = 0; i < offSpringSize; ++i) {
			auto parents = _tournament();
			if (_generateProbability() < CrossoverProbability) {
				_crossover(parents[0], parents[1]);
			}
			for (int j = 0; j < 2; ++j) {
				if (_generateProbability() < MutationProbability) {
					_mutation(parents[j], 1.0 / Population[0].Sequence.size());
				}
				NewPopulation[i * 2 + j] = _createIndividual(parents[j]);
			}
		}
		sort(Population.begin(), Population.end(), [&](const Individual& l, const Individual& r) {
			return _compare(l, r);
		});
		sort(NewPopulation.begin(), NewPopulation.end(), [&](const Individual& l, const Individual& r) {
			return _compare(l, r);
		});
		NewPopulation[PopulationSize - 2] = Population[0];
		NewPopulation[PopulationSize - 1] = Population[1];
		swap(Population, NewPopulation);
	}
	sort(Population.begin(), Population.end(), [&](const Individual& l, const Individual& r) {
		return _compare(l, r);
	});
	BestLeak = Population[0].Fidelity;
	BestAngle = Population[0].OptimizedTheta;
	BestSequence = Population[0].Sequence;
	BestNumberOfCycles = Population[0].NumberOfCycles;
	return omp_get_wtime() - start;
}

void GeneticAlgorithm::reset() {
	BestLeak = numeric_limits<double>::max();
	BestAngle = -1;
	BestSequence.clear();
}

double GeneticAlgorithm::getLeak() {
	return BestLeak;
}

double GeneticAlgorithm::getAngle() {
	return BestAngle;
}

vector<int> GeneticAlgorithm::getSequence() {
	return BestSequence;
}

int GeneticAlgorithm::getNumberOfCycles() {
	return BestNumberOfCycles;
}

void GeneticAlgorithm::_createPopulation(const vector<int>& sequence) {
	Population.clear();
	Population.reserve(2 * sequence.size() + 1);
	Population.push_back(_createIndividual(sequence));
	for (int i = 0; i < sequence.size() * 2; i++) {
		Population.push_back(_createIndividual(Kernel.ChangingOneElement(sequence, i / 2, i % 2)));
	}
	PopulationSize = Population.size();
}

Individual GeneticAlgorithm::_createIndividual(const vector<int>& sequence) {
	vector<int> cur_seq;
	cur_seq.reserve(sequence.size() * Config.NumberOfCycles);
	
	pair<double, double> Q = {INT_MAX, INT_MAX};
	int NumberOfCycles = -1;
	double Theta = INT_MAX;
	
	for (int len = 1; len <= Config.NumberOfCycles; ++len) {
		for (size_t j = 0; j < sequence.size(); ++j) {
			cur_seq.push_back(sequence[j]);
		}
		double OptimizedTheta = _optimizedTheta(cur_seq);
		double F = _fidelity(cur_seq, OptimizedTheta);
		pair<double, double> cur_Q = { fabs(OptimizedTheta - Config.NeededAngle), F };
		if (cur_Q < Q) {
			Q = cur_Q;
			Theta = OptimizedTheta;
			NumberOfCycles = len;
		}
	}
	return Individual(sequence, Theta, Q.second, NumberOfCycles);
}

double GeneticAlgorithm::_optimizedTheta(const vector<int>& Sequence) {
	return Kernel.NewThetaOptimizer(Sequence, Config.Theta);
}

double GeneticAlgorithm::_fidelity(const vector<int>& Sequence, double Theta) {
	return Kernel.Fidelity(Sequence, Theta);
}

vector<vector<int>> GeneticAlgorithm::_tournament() {
	int index1, index2, index3;
	do {
		index1 = _generateInt(0, PopulationSize - 1);
		index2 = _generateInt(0, PopulationSize - 1);
		index3 = _generateInt(0, PopulationSize - 1);
	} while (index1 == index2 || index1 == index3 || index2 == index3);
	vector<Individual> indices = { Population[index1], Population[index2], Population[index3]};
	sort(indices.begin(), indices.end(), [&](const Individual& l, const Individual& r) {
		return _compare(l, r);
	});
	vector<vector<int>> ans = { indices[0].Sequence, indices[1].Sequence };
	return ans;
}

double GeneticAlgorithm::_generateProbability() {
	uniform_int_distribution<> dist(0, numeric_limits<int>::max());
	return 1.0 * dist(RandomGenerator) / numeric_limits<int>::max();
}

int GeneticAlgorithm::_generateInt(int l, int r) {
	uniform_int_distribution<> dist(l, r);
	return dist(RandomGenerator);
}

void GeneticAlgorithm::_crossover(vector<int>& a, vector<int>& b) {
	int index = _generateInt(2, (int)b.size() - 3);
	for (size_t i = index; i < a.size(); ++i) {
		swap(a[i], b[i]);
	}
}

void GeneticAlgorithm::_mutation(vector<int>& a, double p) {
	vector<int> tmp;
	tmp.reserve(3);
	for (int i = 0; i < a.size(); i++) {
		if (_generateProbability() < p) {
			if (Config.Type == 2) tmp = { 0, 1 };
			else tmp = { -1, 0, 1 };
			tmp.erase(find(tmp.begin(), tmp.end(), a[i]));
			if (tmp.size() == 1) a[i] = tmp[0];
			else a[i] = tmp[_generateProbability() * 2 < 1.0];
		}
	}
}

bool GeneticAlgorithm::_compare(const Individual& a, const Individual& b) {
	double left_abs = fabs(a.OptimizedTheta - Config.NeededAngle);
	double right_abs = fabs(b.OptimizedTheta - Config.NeededAngle);
	if (fabs(left_abs - right_abs) < Config.FidelityUpperBound) {
		return a.Fidelity < b.Fidelity;
	}
	return left_abs < right_abs;
}
