#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <complex>
#include <algorithm>
#include <omp.h>
#include <map>
#include "GeneticAlgorithm.h"

using namespace std;

const int TYPE = 3;

int main() {
	setlocale(LC_ALL, "Rus");
	cout.precision(20);
	
	auto program_start = omp_get_wtime();
	const double tstep = 5e-14;
	const double w01 = PI * 8e9;
	const double w12 = w01 - PI * 5e8;
	const double wt = PI * 5e10;
	const double w = 4e-12;
	const double T = PI / wt * 2;
	double Theta = 0.001;
	int CellsNumber = 120;
	const int NumberOfCycles = 1;

	const double Len = (double)CellsNumber * T * NumberOfCycles;
	const double NeededAngle = 0.024;
	const double AngleThreshold = 0.0005;

	Kernel<TYPE> kernel(tstep, w01, w12, wt, w, T);

	auto Amps = kernel.CreateAmpThresholds(CellsNumber);

	double A1 = Amps[0];
	// double A2 = Amps[Amps.size() / 2 - 1];

	cout << "Берём амплитуду А1 = " << A1 << '\n';
	vector<int> SignalString = kernel.CreateStartSCALLOP(CellsNumber, A1);
	double NewTheta = kernel.NewThetaOptimizer(SignalString, Theta);
	cout << "Theta = " << NewTheta << '\n';
	vector<int> A1Sequence;
	double A1F, A1Angle;
	
	int maxIter = 500;
	double crossover_probability = 0.9;
	double mutation_probability = 0.9;

	CalculationDescriptor desc(w01, w12, wt, w, T, Len, tstep, NewTheta, NeededAngle, NumberOfCycles);

	GeneticAlgorithm<TYPE> algo(SignalString, crossover_probability, mutation_probability, maxIter, desc);
	// vector<vector<int>> start_population = {...};
	// GeneticAlgorithm<TYPE> algo(start_population, crossover_probability, mutation_probability, maxIter, desc);
	algo.run();
	A1Sequence = algo.getSequence();
	A1F = algo.getLeak();
	A1Angle = algo.getAngle();
	
	cout << "Для амлитуды А1 = " << A1 << "; угол Theta = " << A1Angle << '\n';
	cout << "Желаемый угол Th0 = " << NeededAngle << '\n';

	cout << "Конечная последовательность\n";
	kernel.WriteSequence(A1Sequence);
	cout << "Утечка F = " << A1F << '\n';
	auto program_end = omp_get_wtime();
	auto program_duration = program_end - program_start;
	cout << "Время работы всей программы = " << program_duration << " сек.\n";
	return 0;
}