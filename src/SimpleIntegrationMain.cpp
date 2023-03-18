#include <algorithm>
#include <cmath>
#include <chrono>
#include <complex>
#include <iostream>
#include <functional>
#include <fstream>
#include <limits.h>
#include <omp.h>
#include <map>
#include <vector>
#include "GeneticAlgorithm.h"
#include <numeric>

using namespace std;

vector<int> convertsSequenceToVector(const string& s, int cycles) {
	vector<int> ans;
	for (int i = 0; i < s.size(); ++i) {
		if (s[i] == '-') {
			ans.push_back(-(s[i + 1] - '0'));
			i++;
		}
		else {
			ans.push_back(s[i] - '0');
		}
	}
	int sz = ans.size();
	for (int i = 0; i < cycles - 1; ++i) {
		for (int j = 0; j < sz; ++j) {
			ans.push_back(ans[j]);
		}
	}
	return ans;
}

int main(int argc, char** argv) {
	double w01Coeff = 4.73047;
	double w01 = w01Coeff * 2 * PI * 1e9;
	const double w12 = w01 - 0.25 * 2 * PI * 1e9;
	const double wt = 25 * 2 * PI * 1e9;
	const double w = 4e-12;
	const double T = PI / wt * 2;
	const double theta = 0.032;
	const double tstep = 1e-15;
	int numberOfCycles = 8;

	Kernel kernel(tstep, w01, w12, wt, w, T);
	string sequence = "1110111100011000100001000110001111011";
	vector<int> seq = convertsSequenceToVector(sequence, numberOfCycles);

	auto res = kernel.Fidelity(seq, theta);
	cout.precision(20);
	cout << "w01 = " << w01Coeff << '\n';
	cout << "Leak = " << res.leak << '\n';
	cout << "F = " << 1 - res.fidelity << '\n';
	return 0;
}