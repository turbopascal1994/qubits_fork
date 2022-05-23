#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <complex>
#include <algorithm>
#include <map>
#include <omp.h>
#include "kernel.h"
#include "linalg.h"
#include "compareDouble.h"

using namespace std;

int main(){
	setlocale(LC_ALL, "Rus");
	cout.precision(20);

	auto program_start = chrono::high_resolution_clock::now();
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

	Kernel<3> kernel(tstep, w01, w12, wt, w, T);

	auto Amps = kernel.CreateAmpThresholds(CellsNumber);

	double A1 = Amps[0];
	double A2 = Amps[Amps.size() / 2 - 1];
	
	cout << "Берём наименьшую амплитуду А1.\n";
	vector<int> SignalString = kernel.CreateStartSCALLOP(CellsNumber, A1);
	double NewTheta = kernel.NewThetaOptimizer(SignalString, Theta);
	cout << NewTheta << '\n';
	cout << "A1 = " << A1 << '\n';

	vector<int> A1Sequence;
	double A1F, A1Angle;
	auto t_start = chrono::high_resolution_clock::now();
	kernel.GenSearch(SignalString, NewTheta, A1Sequence, A1F, A1Angle);
	auto t_end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<std::chrono::seconds>(t_end - t_start).count();
	cout << "Время вычисления алгоритма = " << duration << " сек.\n";
	
	Theta = 0.001;
	cout << "Берём наибольшую амплитуду A2.\n";
	SignalString = kernel.CreateStartSCALLOP(CellsNumber, A2);
	NewTheta = kernel.NewThetaOptimizer(SignalString, Theta);
	cout << "A2 = " << A2 << '\n';
	vector<int> A2Sequence;
	double A2F, A2Angle;
	kernel.GenSearch(SignalString, NewTheta, A2Sequence, A2F, A2Angle);


	// Interpolation
	cout << "Интерполируем амплитуду до нужной, А3.\n";
	auto b_int = (A2Angle - A1Angle) / (A2 - A1);
	auto a_int = A1Angle - b_int * A1;
	auto NewAmpRough = (NeededAngle - a_int) / b_int;

	vector<double> AmpChecks(Amps.size());
	for(size_t i = 0;i < AmpChecks.size();++i){
		AmpChecks[i] = abs(Amps[i] - NewAmpRough);
	}
	int MinAmpCheckN = min_element(AmpChecks.begin(), AmpChecks.end()) - AmpChecks.begin();
	double NewAmp = Amps[MinAmpCheckN];
	Theta = 0.001;

	SignalString = kernel.CreateStartSCALLOP(CellsNumber, NewAmp);
	NewTheta = kernel.NewThetaOptimizer(SignalString, NewTheta);
	vector<int> A3Sequence;
	double A3F, A3Angle;
	kernel.GenSearch(SignalString, NewTheta, A3Sequence, A3F, A3Angle);
	cout << "Для новой амлитуды А3 = " << NewAmp << "; угол Th3 = " << A3Angle << '\n';
	cout << "Желаемый угол Th0 = " << NeededAngle << '\n';
	
	vector<int> NewSequence;
	double NewF = A3F, NewAngle = A3Angle;
	if (Equal(NeededAngle, NewAngle)){
		cout << "Изменение длины последовательности ничего не даст\n";
		NewSequence = A3Sequence;
	}
	else {
		cout << "Теперь меняем длину последовательности М.\n";
		if (Less(NewAngle, NeededAngle + AngleThreshold)){
			while(Not_equal(NewAngle, NeededAngle + AngleThreshold)){
				CellsNumber--;
				cout << "M = " << CellsNumber << '\n';
				SignalString = kernel.CreateStartSCALLOP(CellsNumber, NewAmp);
				kernel.GenSearch(SignalString, NewTheta, NewSequence, NewF, NewAngle);
				if (Equal(NewAngle, NeededAngle) && Less(NewF, 1e-4)){
					break;
				}
			}
		} else if (More(NewAngle, NeededAngle - AngleThreshold)){
			while (Not_equal(NewAngle, NeededAngle - AngleThreshold)) {
				CellsNumber++;
				cout << "M = " << CellsNumber << '\n';
				SignalString = kernel.CreateStartSCALLOP(CellsNumber, NewAmp);
				kernel.GenSearch(SignalString, NewTheta, NewSequence, NewF, NewAngle);
				if (Equal(NewAngle, NeededAngle) && Less(NewF, 1e-4)) {
					break;
				}
			}
		}
	}

	cout << "Желаемый угол Th0 = " << NeededAngle << '\n';
	cout << "Новый угол ThN = " << NewAngle << '\n';
	cout << "Конечная последовательность\n";
	kernel.WriteSequence(NewSequence);
	cout << "Утечка F = " << NewF << '\n';
	auto program_end = chrono::high_resolution_clock::now();
	double program_duration = chrono::duration_cast<std::chrono::milliseconds>(program_end - program_start).count();
	cout << "Время работы программы = " << fixed << setprecision(20) << program_duration / 1000 << " сек.\n";
	return 0;
}