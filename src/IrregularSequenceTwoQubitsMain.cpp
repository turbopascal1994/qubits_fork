#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <fstream>
#include <math.h>
#include <complex>
#include <mkl.h>
#include <omp.h>
#include "linalg.h"
#include <map>
// #include "TwoCubitsKernel.h"

using namespace std;

using TYPE = double;
const TYPE PI = acos(-1.0);

void fillIdentity(vector<complex<TYPE>>& A, int dim) {
	for (size_t i = 0; i < dim; ++i) {
		A[i * dim + i] = 1;
	}
}

void matadd(const vector<complex<TYPE>>& A, const vector<complex<TYPE>>& B, vector<complex<TYPE>>& Res) {
	for (size_t i = 0; i < A.size(); ++i) {
		Res[i] = A[i] + B[i];
	}
}

void matsub(const vector<complex<TYPE>>& A, const vector<complex<TYPE>>& B, vector<complex<TYPE>>& Res) {
	for (size_t i = 0; i < A.size(); ++i) {
		Res[i] = A[i] - B[i];
	}
}

void vsMul(const vector<complex<TYPE>>& A, complex<TYPE> b, vector<complex<TYPE>>& Res) {
	for (size_t i = 0; i < A.size(); ++i) {
		Res[i] = A[i] * b;
	}
}

void mvMul(const vector<complex<TYPE>>& A, const vector<complex<TYPE>>& B, vector<complex<TYPE>>& Res) {
	for (size_t i = 0; i < Res.size(); ++i) {
		Res[i] = 0;
		for (size_t j = 0; j < B.size(); ++j) {
			Res[i] += A[i * B.size() + j] * B[j];
		}
	}
}

void kMul(const vector<complex<TYPE>>& A, const vector<complex<TYPE>>& B, vector<complex<TYPE>>& Res, int dim) {
	for (size_t i = 0; i < dim * dim; ++i) {
		for (size_t j = 0; j < dim * dim; ++j) {
			// Res[i][j] = A[i / Size][j / Size] * B[i % Size][j % Size];
			Res[i * dim * dim + j] = A[(i / dim) * dim + j / dim] * B[(i % dim) * dim + j % dim];
		}
	}
}

complex<TYPE> trace(const vector<complex<TYPE>>& A, int dim) {
	complex<TYPE> trace = 0;
	for (int i = 0; i < dim; ++i) {
		trace += A[i * dim + i];
	}
	return trace;
}

void ctranspose(const vector<complex<TYPE>>& A, vector<complex<TYPE>>& Res, int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
			Res[i * dim + j] = conj(A[j * dim + i]);
			Res[j * dim + i] = conj(A[i * dim + j]);
		}
	}
}

vector<int> OperatorGridIrregular(
	TYPE w,
	TYPE T1,
	TYPE T2,
	TYPE phi,
	int waitq1,
	int waitq2,
	string str1,
	string str2,
	TYPE step
) {
	/*
	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % combines two grids and decides which combined operator  %
    % should be applied on each grid step                     %
    % 1 grid step = 1 grs                                     %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
	*/
	int Nw = ceil(w / step); // pulse width(in grs)
	int NT1 = ceil(T1 / step) - Nw; // distance b / w pulses on Q1(in grs)
	int NT2 = ceil(T2 / step) - Nw; // distance b / w pulses on Q2(in grs)
	vector<int> sz(Nw, 0); // array for 0 pulse
	vector<int> sp(Nw, 1); // array for + 1 pulse
	vector<int> sm(Nw, -1); // array for - 1 pulse
	vector<int> s0_1(NT1, 0); // array for distance b / w pulses on Q1
	vector<int> s0_2(NT2, 0); // array for distance b / w pulses on Q2
	vector<int> arr1, arr2;
	
	// phase between pulses(in grs)
	/*if (phi > 0)
		phi_arr = zeros(1, phi);
	end*/

	// placing pulsesand the distance b / w them
	for (int j = 0; j < str1.size(); ++j) {
		if (str1[j] == '0') {
			arr1.insert(arr1.end(), sz.begin(), sz.end());
		}
		else if (str1[j] == '1') {
			arr1.insert(arr1.end(), sp.begin(), sp.end());
		}
		else {
			j += 1;
			arr1.insert(arr1.end(), sm.begin(), sm.end());
		}
		if (j + 1 != str1.size()) {
			arr1.insert(arr1.end(), s0_1.begin(), s0_1.end());
		}
	}
	for (int j = 0; j < str2.size(); ++j) {
		if (str2[j] == '0') {
			arr2.insert(arr2.end(), sz.begin(), sz.end());
		}
		else if (str2[j] == '1') {
			arr2.insert(arr2.end(), sp.begin(), sp.end());
		}
		else {
			j += 1;
			arr2.insert(arr2.end(), sm.begin(), sm.end());
		}
		if (j + 1 != str2.size()) {
			arr2.insert(arr2.end(), s0_2.begin(), s0_2.end());
		}
	}
	
	// wait time after pulses(in grs)
	for (int i = 0; i < waitq1; i++) arr1.push_back(0);
	for (int i = 0; i < waitq2; i++) arr2.push_back(0);
	
	// equalizing pulse strings by adding zeros to the lesser one
	while (arr1.size() > arr2.size()) {
		arr2.push_back(0);
	}
	while (arr1.size() < arr2.size()) {
		arr1.push_back(0);
	}
	vector<int> arr;
	/*
	combining into one string
    legend:
			0 = -1-1
			1 = 0-1
			2 = 1-1
			3 = -10
			4 = 00
			5 = 10
			6 = -11
			7 = 01
			8 = 11
	*/
	for (int i = 0; i < arr1.size(); i++) {
		arr.push_back((arr1[i] + 1) + (arr2[i] + 1) * 3);
	}
	return arr;
}

map<string, TYPE> SimulateIrregular(
	int N,
	TYPE w1, 
	TYPE w2, 
	TYPE mu1, 
	TYPE mu2, 
	TYPE g, 
	TYPE Cq1, 
	TYPE Cq2, 
	TYPE Cc1, 
	TYPE Cc2,
	TYPE wg1, 
	TYPE wg2, 
	TYPE tau, 
	TYPE phi, 
	int waitq1, 
	int waitq2,
	string str1, 
	string str2, 
	TYPE tstep, 
	string init, 
	string operation
) {
	const TYPE h = 1.054e-34;// planck constant
	const TYPE F0 = 2.06e-15;// magnetic flux quantum
	const complex<TYPE> zero = { 0, 0 };

	vector<complex<TYPE>> a1(N * N), 
		a2(N * N), 
		Identity(N * N), 
		lIdentity(N * N * N * N),
		aa, V1(N * N), V2(N * N),
		HQ1(N * N), HQ2(N * N);
	fill(a1.begin(), a1.end(), zero);
	fill(a2.begin(), a2.end(), zero);
	// second order quantizaion
	for (int i = 1; i < N; i++) {
		a1[i * N + i - 1] = { sqrt(i), 0 };
		a2[(i - 1) * N + i] = { sqrt(i), 0 };
	}
	fillIdentity(Identity, N);
	fillIdentity(lIdentity, N * N);
	aa = linalg::matmul(a1, a2, N, N, N, N, N, N);
	
	// field operator
	TYPE V0 = F0 / tau; // Voltage
	TYPE Amp1 = Cc1 * V0 * sqrt(h * w1 / (2 * Cq1)); //first generator amplitude
	TYPE Amp2 = Cc2 * V0 * sqrt(h * w2 / (2 * Cq2)); //second generator amplitude

	vector<complex<TYPE>> tmp1(N * N),
		tmp2(N * N), tmp3(N * N), tmp4(N * N), tmp5(N * N);
	matsub(a2, a1, tmp1);
	vsMul(tmp1, complex<TYPE>{0, Amp1}, V1);
	vsMul(tmp1, complex<TYPE>{0, Amp2}, V2);


	// now to the hamiltonians
	int L = N * N;
	vector<complex<TYPE>> 
		ltmp1(L * L), ltmp2(L * L), 
		ltmp3(L * L), ltmp4(L * L), ltmp5(L * L),
		H00(L * L);
	vsMul(aa, h * w1, tmp1);
	vsMul(aa, h * w2, tmp3);

	matsub(aa, Identity, tmp5);
	tmp4 = linalg::matmul(aa, tmp5, N, N, N, N, N, N);
	vsMul(tmp4, h * mu1 / 2, tmp2);
	matsub(tmp1, tmp2, HQ1);
	vsMul(tmp4, h * mu2 / 2, tmp2);
	matsub(tmp3, tmp2, HQ2);

	matadd(a1, a2, tmp2);
	kMul(tmp2, tmp2, ltmp1, N);
	vector<complex<TYPE>> Hint(L * L);
	vsMul(ltmp1, h * g, Hint);
	

	// no field hamiltonian
	kMul(HQ1, Identity, ltmp2, N);
	kMul(Identity, HQ2, ltmp3, N);
	matadd(ltmp2, ltmp3, ltmp1);
	matadd(ltmp1, Hint, H00);

	
	// eigens
	vector<complex<TYPE>> 
		EigVectorsL(L * L),
		EigVectorsR(L * L),
		EigValues(L * L);
	linalg::eig(H00, EigVectorsL, EigVectorsR, EigValues, L);
	
	vector<int> IndexEigValuesAndVectors(N * N);
	iota(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), 0);
	sort(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), [&](int el1, int el2) {
		return EigValues[el1].real() < EigValues[el2].real();
	});
	auto getEigVector = [&](int index) {
		vector<complex<TYPE>> vec(L);
		for (int i = 0; i < L; ++i) {
			vec[i] = EigVectorsR[i * L + IndexEigValuesAndVectors[index]];
		}
		return vec;
	};
	vector<complex<TYPE>> WF00, WF10, WF01, WF11, WF20, WF02;

	WF00 = getEigVector(0);
	WF10 = getEigVector(1);
	WF01 = getEigVector(2);
	if (N == 2) {
		WF11 = getEigVector(3);
	}
	else if (N == 3) {
		WF20 = getEigVector(3);
		WF02 = getEigVector(4);
		WF11 = getEigVector(5);
	}

	vector<complex<TYPE>>
		H10(L * L),
		H01(L * L),
		H11(L * L),
		Hm10(L * L),
		H0m1(L * L),
		Hm11(L * L),
		H1m1(L * L),
		Hm1m1(L * L);

	kMul(V1, Identity, ltmp1, N);
	kMul(Identity, V2, ltmp2, N);
	vsMul(V1, -1, tmp1);
	kMul(tmp1, Identity, ltmp3, N);
	vsMul(V2, -1, tmp2);
	kMul(Identity, tmp2, ltmp4, N);

	matadd(H00, ltmp1, H10);
	matadd(H00, ltmp2, H01);
	matadd(H10, ltmp2, H11);
	matadd(H00, ltmp3, Hm10);
	matadd(H00, ltmp4, H0m1);
	matadd(Hm10, ltmp2, Hm11);
	matadd(H10, ltmp4, H1m1);
	matadd(Hm10, ltmp4, Hm1m1);
	

	vector<vector<complex<TYPE>>> operators = {
		linalg::getUMatrix(lIdentity, Hm1m1, tstep, h, L),
		linalg::getUMatrix(lIdentity, H0m1, tstep, h, L),
		linalg::getUMatrix(lIdentity, H1m1, tstep, h, L),
		linalg::getUMatrix(lIdentity, Hm10, tstep, h, L),
		linalg::getUMatrix(lIdentity, H00, tstep, h, L),
		linalg::getUMatrix(lIdentity, H10, tstep, h, L),
		linalg::getUMatrix(lIdentity, Hm11, tstep, h, L),
		linalg::getUMatrix(lIdentity, H01, tstep, h, L),
		linalg::getUMatrix(lIdentity, H11, tstep, h, L)
	};
	
	vector<int> PulseString = OperatorGridIrregular(
		tau, 2 * PI / wg1, 2 * PI / wg2, phi, waitq1, waitq2, str1, str2, tstep
	);
	vector<complex<TYPE>> U(L * L);
	fillIdentity(U, L);
	for (int i = 0; i < PulseString.size(); i++) {
		auto& UPulse = operators[PulseString[i]];
		U = linalg::matmul(UPulse, U, L, L, L, L, L, L);
	}
	
	vector<complex<TYPE>> WF;
	if (init == "00") {
		WF = WF00;
	}
	else if (init == "01") {
		WF = WF01;
	}
	else if (init == "10") {
		WF = WF10;
	}
	else if (init == "11") {
		WF = WF11;
	}
	else if (init == "20") {
		WF = WF20;
	}
	else if (init == "02") {
		WF = WF02;
	}
	mvMul(U, WF, tmp1);
	swap(WF, tmp1);

	auto getProbability = [&](const vector<complex<TYPE>>& eigVector) {
		complex<TYPE> dot_product = 0;
		for (int j = 0; j < L; ++j) {
			dot_product += conj(eigVector[j]) * WF[j];
		}
		return norm(dot_product);
	};

	map<string, TYPE> probs = {
		{"00", getProbability(WF00)},
		{"10", getProbability(WF10)},
		{"01", getProbability(WF01)},
		{"11", getProbability(WF11)}
	};
	if (N == 3) {
		probs["20"] = getProbability(WF20);
		probs["02"] = getProbability(WF02);
	}
	

	// THIS IS THE PART THAT CALCULATES FIDELITY(MAY BE WRONG)
	// 1Q ideal gate matrices
	TYPE dth = PI / 2;
	vector<complex<TYPE>> Ypi2 = {
		cos(dth / 2), -sin(dth / 2), 0,
		sin(dth / 2), cos(dth / 2), 0,
		0, 0, 1
	};
	vector<complex<TYPE>> Ypi = linalg::matmul(Ypi2, Ypi2, N, N, N, N, N, N);
	// 2Q ideal gate matrices
	vector<complex<TYPE>>
		Y00(L * L), Yh0(L * L), Y0h(L * L),
		Y10(L * L), Y01(L * L), Yhh(L * L),
		Y1h(L * L), Yh1(L * L), Y11(L * L);
	kMul(Identity, Identity, Y00, N);
	kMul(Ypi2, Identity, Yh0, N);
	kMul(Identity, Ypi2, Y0h, N);
	kMul(Ypi, Identity, Y10, N);
	kMul(Identity, Ypi, Y01, N);
	kMul(Ypi2, Ypi2, Yhh, N);
	kMul(Ypi, Ypi2, Y1h, N);
	kMul(Ypi2, Ypi, Yh1, N);
	kMul(Ypi, Ypi, Y11, N);

	vector<complex<TYPE>> Uid;
	if (operation == "00") Uid = Y00;
	else if (operation == "01") Uid = Y01;
	else if (operation == "10") Uid = Y10;
	else if (operation == "11") Uid = Y11;
	else if (operation == "h0") Uid = Yh0;
	else if (operation == "0h") Uid = Y0h;
	else if (operation == "h1") Uid = Yh1;
	else if (operation == "1h") Uid = Y1h;
	else if (operation == "hh") Uid = Yhh;

	vector<complex<TYPE>> conj_Uid(Uid.size());
	ctranspose(Uid, conj_Uid, L);
	auto M = linalg::matmul(conj_Uid, U, L, L, L, L, L, L);
	vector<complex<TYPE>> conj_M(Uid.size());
	ctranspose(M, conj_M, L);
	auto MMconj = linalg::matmul(M, conj_M, L, L, L, L, L, L);
	TYPE F = (abs(trace(MMconj, L)) + norm(trace(M, L))) / (N * N * (N * N + 1));
	probs["F"] = F;
	return probs;
}

int main() {
	double time = omp_get_wtime();
	int N = 3; // кол-во уровней кубита
	TYPE val = 2 * PI * 1e9;
	TYPE tstep = 5e-14; // time grid step
	// main qubit frequencies
	TYPE w1 = 5.0 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	TYPE w2 = 5.2 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	// anharmonicities
	TYPE mu1 = 0.25 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	TYPE mu2 = 0.4 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	TYPE g = 0.02 * (2 * PI) * 1e9; // параметр взаимодействия между кубитами
	
	// qubit capacities
	TYPE Cq1 = 1e-12;
	TYPE Cq2 = 1e-12;

	// connection capacities
	TYPE Cc1 = 4.9e-16;
	TYPE Cc2 = 4e-16;

	// pulse generation frequencies
	TYPE wg1 = 25 * (2 * PI) * 1e9;
	TYPE wg2 = 25 * (2 * PI) * 1e9;
	TYPE tau = 4 * 1e-12; // Длительность импульса
	TYPE phi = 0; // phase (number of grid steps paused on Q2)

	// wait time after pulse
	int waitq1 = 0;
	int waitq2 = 0;

	string str1 = "11-1-1-1110-1011-1-1-111-1-1-111-1-1110-1-1110-1-1110-1-111-10-1111-1-1110-1-1111-1-1111-1-1111-1-1011-10111-1-1011-1-1-111-1-1-111-10110-1-111-1-1-101-1-1-111-1-1-1110-1-11";
	string str2;

	string init = "00"; // initial condition
	string operation = "h0"; // required operation (for fidelity calculation)

	auto res = SimulateIrregular(
		N, w1, w2, mu1, mu2, g, Cq1, Cq2, Cc1, Cc2,
		wg1, wg2, tau, phi, waitq1, waitq2,
		str1, str2, tstep, init, operation
	);
	for (auto& i : res) {
		cout << fixed << setprecision(20) << i.first << ' ' << i.second << '\n';
	}
	return 0;
}



