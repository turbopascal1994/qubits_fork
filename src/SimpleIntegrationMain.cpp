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
#include "ConstantsDescriptor.h"
#include "kernel.h"
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
	double w01Coeff = 5;
	double w01 = w01Coeff * 2 * PI * 1e9;
	const double w12 = w01 - 0.25 * 2 * PI * 1e9;
	const double wt = 35 * 2 * PI * 1e9;
	const double T = PI / wt * 2;
	const double theta = 0.024;
	const double tstep = 5e-14;
	const double w = tstep;
	int numberOfCycles = 4;

	ConstantsDescriptor config(w01, w12, wt, w, T, tstep, theta, theta, numberOfCycles, 0.0001, 3);
	Kernel kernel(config);
	string sequence = "000000-11000010000-1-1110000-1000-1-1-111110-1-111100-111110110-1110000-1100-1-1-1000111110011100-1111111111-1111110-111-1-111000-1-11111110111000-1111111-1110-1-111111111100000-1-1-11110-1-1-11000000-11100-1-1";
	vector<int> seq = convertsSequenceToVector(sequence, numberOfCycles);
	cout << seq.size() << '\n';

	ofstream out;
	out.open("our_solve.txt");
	auto res = kernel.Fidelity(seq, config.neededAngle, &out);
	cout.precision(20);
	cout << "w01 = " << w01Coeff << '\n';
	cout << "Leak = " << res.leak << '\n';
	cout << "F = " << 1 - res.fidelity << '\n';
	return 0;
}