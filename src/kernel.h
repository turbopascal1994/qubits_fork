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

	struct IntegrationResult {
		double level0, level1, level2;
		vector<complex<double>> state;
	};

	double tstep, w01, w12, wt, w, T;
	vector<complex<double>> Id, H0, EigVec, EigVecRight, EigVal, WF1, WF2, WF3, Hmatrix, InitStates, IdealStates;

	void precalcFidelityAndRotateCheck(double tstep, double w01, double w12);
	IntegrationResult integration(
		int N,
		const vector<int>& InputSequence,
		const vector<complex<double>>& UTZero,
		const vector<complex<double>>& UTMinus,
		const vector<complex<double>>& UTPlus,
		vector<complex<double>>& WF
	);
	double RotateCheck(const vector<int>& InputSequence, double Theta);
	void prepareUMatrices(double Theta, vector<complex<double>>& UTZero, vector<complex<double>>& UTMinus, vector<complex<double>>& UTPlus);
public:

	struct FidelityResult {
		double fidelity, leak;
	};

	Kernel(
		double tstep,
		double w01,
		double w12,
		double wt,
		double w,
		double T): tstep(tstep), w01(w01), w12(w12), wt(wt), w(w), T(T) {
		precalcFidelityAndRotateCheck(tstep, w01, w12);
	}

	vector<int> CreateStartSCALLOP(int M, double AmpThreshold, int TYPE);

	vector<double> CreateAmpThresholds(const int M);

	FidelityResult Fidelity(const vector<int>& SignalString, double Theta);

	double NewThetaOptimizer(const vector<int>& InputSequence, double UnOptTheta);

	void WriteSequence(ostream& fout, const vector<int>& InputString, char End = '\n');

	vector<int> ChangingOneElement(vector<int> InputSequence, int index, int type, int TYPE);
};
