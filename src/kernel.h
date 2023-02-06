#pragma once
#include <complex>
#include <tuple>
#include <omp.h>
#include "linalg.h"

using namespace std;
static const double PI = acos(-1.0);
static const double h = 1.054e-34;
static const double C1 = 1e-12;
static const double F0 = 2.06e-15;

class Kernel {
	int TYPE;
	double tstep, w01, w12, wt, w, T;
	vector<complex<double>> Id, H0, EigVec, EigVal, WF1, WF2, WF3, Hmatrix, InitStates;

	void precalcFidelityAndRotateCheck(double tstep, double w01, double w12);
	double calcProbability(
		int N,
		const vector<int>& InputSequence,
		const int CyclePlusMinusSteps,
		const int CycleZeroSteps,
		const vector<complex<double>>& UPlus,
		const vector<complex<double>>& UPlus2,
		const vector<complex<double>>& UMinus,
		const vector<complex<double>>& UMinus2,
		const vector<complex<double>>& UZero,
		vector<complex<double>>& WF,
		const vector<complex<double>>& WF1);
	double RotateCheck(const vector<int>& InputSequence, double Theta);
public:
	Kernel(
		double tstep,
		double w01,
		double w12,
		double wt,
		double w,
		double T,
		int Type): tstep(tstep), w01(w01), w12(w12), wt(wt), w(w), T(T), TYPE(Type) {
		precalcFidelityAndRotateCheck(tstep, w01, w12);
	}

	vector<int> CreateStartSCALLOP(int M, double AmpThreshold);

	vector<double> CreateAmpThresholds(const int M);

	tuple<double, double, double, double> Fidelity(const vector<int>& SignalString, double Theta);

	double NewThetaOptimizer(const vector<int>& InputSequence, double UnOptTheta);

	void WriteSequence(ostream& fout, const vector<int>& InputString, char End = '\n');

	vector<int> ChangingOneElement(vector<int> InputSequence, int index, int type);

	void GenSearch(
		vector<int> InputString,
		double StartTheta, vector<int>& BestSequenceOverall,
		double& BestLeakOverall,
		double& BestAngleOverall
	);
};
