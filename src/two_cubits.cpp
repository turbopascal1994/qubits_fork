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

using namespace std;

using TYPE = double;

const TYPE PI = acos(-1.0);

//----------------- Constants ---------------------------- -
//TYPE, parameter::pi = 3.1415926535897932384626433832795
//TYPE, parameter::hh = 1.054 * 1E-34
//TYPE, parameter::F0 = 2.06 * 1E-15
//TYPE, parameter::C1 = 1. * 1E-12
TYPE hh = 1.054 * 1E-34; // приведённая постоянная Планка
TYPE F0 = 2.06 * 1E-15; // квант потока
TYPE C1 = 1.0 * 1E-12; // емкость кубита
//
//
//!integer, parameter::L = 5.
//complex, parameter::Ic = (0., 1.)
complex<TYPE> Ic(0.0, 1.0);
//
//integer, parameter::Nm = 3.
const int Nm = 3; // Число учитываемых уровней одного кубита
//integer, parameter::L = Nm * Nm
const int L = Nm * Nm; // Размерность системы в целом
//TYPE, parameter::w0q1 = 5.12 * (2 * pi) * 1E9 !5.3463
TYPE w0q1 = 5.12 * (2 * PI) * 1e9; // Собственная частота первого кубита
//TYPE, parameter::w0q2 = 5.350 * (2 * pi) * 1E9  !5.1167
TYPE w0q2 = 5.350 * (2 * PI) * 1e9; // Собственная частота второго кубита
//
//TYPE, parameter::wr = 7 * (2 * pi) * 1E9
//TYPE, parameter::g1 = 0.07 * (2 * pi) * 1E9
//TYPE, parameter::g2 = 0.07 * (2 * pi) * 1E9
//TYPE, parameter::g = 0.02 * (2 * pi) * 1E9!Abs((g1 * g2 * (w0q1 + w0q2 - 2 * wr)) / (2 * (w0q1 - wr) * (w0q2 - wr))) !!!0.01 * (2 * pi) * 1E9!0.0017 * (2 * pi) * 1E9!
TYPE g = 0.02 * (2 * PI) * 1e9; // параметр взаимодействия между кубитами
//
//TYPE, parameter::dw = w0q1 - w0q1  !5.1167
//TYPE, parameter::wt = 5.1258 * (2 * pi) * 1E9 !g = 0.015 * (2 * pi) * 1E9
TYPE wt = 5.1258 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
//
//TYPE, parameter::mu1 = -0.353 * (2 * pi) * 1E9
//TYPE, parameter::mu2 = -0.35 * (2 * pi) * 1E9
TYPE mu1 = -0.353 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
TYPE mu2 = -0.35 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита

//integer, parameter::M = 20000
const int M = 20000; // Максимальное число периодов внешнего управляющего поля
//TYPE, parameter::d = 2 * pi / (wt)
TYPE d = 2 * PI / (wt); // Период внешнего управляющего поля
//TYPE, parameter::tau = 4 * 1E-12
TYPE tau = 4 * 1e-12; // Длительность импульса
//
//TYPE, parameter::dr = 4 * 1E-13
//TYPE, parameter::dt = 2 * 1E-13, ddt = dt / 2.
TYPE dt = 2 * 1e-13, ddt = dt / 2.0;
//TYPE, parameter::tmax = 270 * 1E-9
TYPE TMax = 270 * 1e-9; // Время интегрирования
//integer, parameter::Nsteps = tmax / dt
int StepsNumber = (int)(TMax / dt); // Число шагов интегрирования
//
//TYPE, parameter::tetta = 0.0018
TYPE Theta = 0.0018; // Измеренный угол поворота
//TYPE, parameter::V = F0 / tau
TYPE V = F0 / tau; // "Высота" импульса
//TYPE, parameter::tcv = tau / 2.0
//
//!TYPE, parameter::tx = 30 * 1E-9
//!TYPE, parameter::Trise = 50 * 1E-9
//!integer, parameter::irise = Trise / dr
//!TYPE, parameter::sigmax = Trise / 4.0
//!TYPE, parameter::Ax = 0.05 * (2 * pi) * 1E9 !0.00221 Pi / 2 (k = 12)
//
//!TYPE, parameter::tcr = 330 * 1E-9 !128.193 * 1E-9 !41.865 * 1E-9!128.193 * 1E-9
//!TYPE, parameter::tc = 2 * Trise
//!integer, parameter::itcr = tcr / dr
//!!TYPE, parameter::Trise = 15 * 1E-9
//!TYPE, parameter::sigmacr = Trise / 3.0
//!TYPE, parameter::Acr = 0.036 * (2 * pi) * 1E9 !0.0940 Pi / 2 (k = 12)
//!TYPE, parameter::Acon = -0.001944 * (2 * pi) * 1E9 !0.00162 Pi / 2 (k = 12)
//
//integer, parameter::LWMAX = 1000
//
//complex, dimension(Nm, Nm) ::h1, h2, nn, h1h2
vector<complex<TYPE>> H1(Nm* Nm);
vector<complex<TYPE>> H2(Nm* Nm);
vector<complex<TYPE>> H1H2(Nm* Nm);
vector<complex<TYPE>> NN(Nm* Nm);
//complex, dimension(L, L) ::h12, Hint, H01, H02, H0, V1, V2, Hr, VL, VR, W, Rr, Rl, Edd, Ud, dU
// Гамильтонианы для общей системы
vector<complex<TYPE>> H12(L* L); // служебный
vector<complex<TYPE>> HInteration(L* L); // Гамильтониан взаимодействия
vector<complex<TYPE>> H01(L* L); // Гамильтониан первого кубита в контексте общей системы
vector<complex<TYPE>> H02(L* L); // Гамильтониан второго кубита в контексте общей системы
vector<complex<TYPE>> H0(L* L); // Cтационарный гамильтониан двухкубитной системы
vector<complex<TYPE>> EigVectorsL(L* L); // Левые собственные вектора гамильтониана двухкубитной системы
vector<complex<TYPE>> EigVectorsR(L* L); // Правые собственные вектора гамильтониана двухкубитной системы
vector<complex<TYPE>> EigValues(L); // Собственные числа гамильтониана двухкубитной системы
vector<int> IndexEigValuesAndVectors(L); // Перестановка номеров собственных чисел и векторов для упоредочивания их по возрастанию действительной части собственных чисел
vector<complex<TYPE>> V1(L* L); // Сумма операторов рождения и уничтожения для первого кубита в контексте общей системы
vector<complex<TYPE>> V2(L* L); // Сумма операторов рождения и уничтожения для второго кубита в контексте общей системы
vector<complex<TYPE>> Hr(L* L); // Гамильтониан двухкубитной системы на каждом шаге по времени
vector<complex<TYPE>> Rr(L* L);
vector<complex<TYPE>> Rl(L* L);
vector<complex<TYPE>> dU(L* L); // Оператор системы на шаге dt
//complex, dimension(Nm, Nm) ::Ed, H0q1, H0q2
vector<complex<TYPE>> Identity(Nm* Nm);
vector<complex<TYPE>> nIdentity(Nm* Nm);
vector<complex<TYPE>> lIdentity(L* L);
vector<complex<TYPE>> H0q1(Nm* Nm); // Независимые гамильтонианы первого и второго кубитов
vector<complex<TYPE>> H0q2(Nm* Nm);
//TYPE::buf, a11, a12, a13, a14, a15, a16, a17, a18, a19, d2
TYPE d2; // ??? Половина периода ...
//TYPE::Cc, A, Ka, At
TYPE Cc; //
TYPE A; //
TYPE At; //
TYPE Ka; //
//complex, dimension(L, L) ::evea, eveb
//complex, dimension(L) ::eign, fi
int CurEigVectorNumber;
vector<complex<TYPE>> CurEigVector(L); // Собственный вектор гамильтониана двухкубитной системы, для которого выполняется интегрирование
vector<complex<TYPE>> UpdatedVector(L);
//integer, dimension(L) ::ipiv
//integer, dimension(M + 1) ::Am
//!TYPE, dimension(M* Nc + 1) ::Ap
//TYPE, dimension(M + 1) ::Ap
vector<TYPE> Ap(M); // ??? Амплитуды ...
//integer::i, jt, j, s, jj
//TYPE::it, t
TYPE TCurrent;
//integer::info, LWORK
//TYPE             RWORK(2 * L)
//COMPLEX          WORK(LWMAX)
//TYPE, DIMENSION(L) ::wor

vector<complex<TYPE>> tmp1(Nm* Nm);
vector<complex<TYPE>> tmp2(Nm* Nm);
vector<complex<TYPE>> tmp3(Nm* Nm);
vector<complex<TYPE>> tmp4(Nm* Nm);

vector<complex<TYPE>> ltmp1(L* L);
vector<complex<TYPE>> ltmp2(L* L);
vector<complex<TYPE>> ltmp3(L* L);

void FillIdentity(vector<complex<TYPE>>& A, int dim) {
	for (size_t i = 0; i < dim; ++i) {
		A[i * dim + i] = 1;
	}
}

void matadd(const vector<complex<TYPE>>& A, const vector<complex<TYPE>>& B, vector<complex<TYPE>>& Res, int N) {
	/*Res = B;
	vector<complex<TYPE>> Id(N * N);
	FillIdentity(Id, N);
	MKL_Complex16 one = { 1, 0 };
	MKL_Complex16 zero = { 0, 0 };
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &one, A.data(), N, Id.data(), N, &one, Res.data(), N);*/
	for (size_t i = 0; i < A.size(); ++i) {
		Res[i] = A[i] + B[i];
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

void kMul(const vector<complex<TYPE>>& A, const vector<complex<TYPE>>& B, vector<complex<TYPE>>& Res) {
	int dim = sqrt(sqrt(Res.size()));
	for (size_t i = 0; i < dim * dim; ++i) {
		for (size_t j = 0; j < dim * dim; ++j) {
			//      Res[i][j] = A[i / Size][j / Size] * B[i % Size][j % Size];
			Res[i * dim * dim + j] = A[(i / dim) * dim + j / dim] * B[(i % dim) * dim + j % dim];
		}
	}
}

int main() {
	double time = omp_get_wtime();
	Cc = Theta / (sqrtl(2 * wt * F0 * F0) / (sqrtl(hh) * sqrtl(C1)));
	A = 2 * Cc * V * sqrtl((hh * wt) / (C1 * 2));
	d2 = d / 2.0;

	// операторы рождения, уничтожения и их сумма, оператор числа частиц
	complex<TYPE> zero = { 0, 0 };
	fill(H1.begin(), H1.end(), zero);
	fill(H2.begin(), H2.end(), zero);
	fill(H1H2.begin(), H1H2.end(), zero);
	fill(NN.begin(), NN.end(), zero);

	H1[1] = H2[Nm] = 1;
	for (int i = 1; i < Nm - 1; i++) {
		H1[i * Nm + i + 1] = { sqrtl(i + 1), 0 };
		H2[(i + 1) * Nm + i] = { sqrtl(i + 1), 0 };
	}

	NN = linalg::matmul(H2, H1, Nm, Nm, Nm, Nm, Nm, Nm);
	matadd(H1, H2, H1H2, Nm);

	FillIdentity(Identity, Nm);
	FillIdentity(lIdentity, L);
	vsMul(Identity, -1.0, nIdentity);

	// гамильтонианы отдельных кубитов размерностью NxN
	vsMul(NN, complex<TYPE>(hh * w0q1), tmp1);
	matadd(NN, nIdentity, tmp2, Nm);
	tmp3 = linalg::matmul(NN, tmp2, Nm, Nm, Nm, Nm, Nm, Nm);
	vsMul(tmp3, complex<TYPE>(0.5 * mu1 * hh), tmp4);
	matadd(tmp1, tmp4, H0q1, Nm);

	vsMul(NN, complex<TYPE>(hh * w0q2), tmp1);
	matadd(NN, nIdentity, tmp2, Nm);
	tmp3 = linalg::matmul(NN, tmp2, Nm, Nm, Nm, Nm, Nm, Nm);
	vsMul(tmp3, complex<TYPE>(0.5 * mu2 * hh), tmp4);
	matadd(tmp1, tmp4, H0q2, Nm);

	/*linalg::print_matrix("H0q1", Nm, Nm, H0q1, Nm);
	linalg::print_matrix("H0q2", Nm, Nm, H0q2, Nm);*/

	// // Гамильтониан взаимодействия
	kMul(H1H2, H1H2, H12);
	vsMul(H12, complex<TYPE>(hh * g), HInteration);

	// // Гамильтонианы кубитов в контексте общей системы
	kMul(H0q1, Identity, H01);
	kMul(Identity, H0q2, H02);

	// // Cтационарный гамильтониан двухкубитной системы
	matadd(H01, H02, ltmp1, L);
	matadd(ltmp1, HInteration, H0, L);

	//linalg::print_matrix("H0", L, L, H0, L);

	// // Вычисление собственных чисел и векторов стационарного гамильтониана двухкубитной системы
	linalg::eig(H0, EigVectorsL, EigVectorsR, EigValues, L);

	/*linalg::print_matrix("EigVectorsR 1", L, 1, EigVectorsR, L);*/
	/*
	int index = 2;
	cout << EigValues[index] << '\n';
	for (int i = 0; i < L; i++) {
		cout << EigVectorsR[i * L + index] << '\n';
	}
	linalg::print_matrix("EigValues", L, 1, EigValues, 1);*/

	/*for (int index = 0; index < L; ++index) {
		cout << "EigVal = " << EigValues[index] << ";\t";
		TYPE ma = 0;
		for (int i = 0; i < L; ++i) {
			complex<TYPE> productL = 0, productR = EigValues[index] * EigVectorsR[i * L + index];

			for (int j = 0; j < L; ++j) {
				productL += conj(H0[i * L + j]) * EigVectorsR[j * L + index];
			}

			ma = max(ma, (productL - productR).real());
		}
		cout << "max diff = " << ma << '\n';
	}

	for (int index = 0; index < L; ++index) {
		auto cp = H0;
		for (int i = 0; i < L; ++i) {
			cp[i * L + i] -= EigValues[index];
		}
		cout << "EigVal = " << EigValues[index] << "\n";
		vector<complex<TYPE>> b(L, 0);
		auto x = linalg::solve_equations(cp, b, L);
		for (auto& i : x) cout << i << ' ';
		cout << '\n';
	}
	return 0;*/

	// // Сортировка собственных чисел по возрастанию действительной части
	iota(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), 0);

	sort(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), [&](int el1, int el2) {
		return EigValues[el1].real() < EigValues[el2].real();
	});

	/*cout << "EVEA\n";
	for(int j = 0;j < L;j++, cout << '\n')
	for (int i = 0; i < L; i++) {
		cout << EigVectorsR[j * L + IndexEigValuesAndVectors[i]].real() << '\t';
	}
	cout << '\n';*/
	

	// kMul(H1H2, Identity, V1);
	kMul(Identity, H1H2, V2);

	CurEigVectorNumber = 2;
	for (int i = 0; i < L; i++) {
		CurEigVector[i] = EigVectorsR[i * L + IndexEigValuesAndVectors[CurEigVectorNumber]];
	}

	// linalg::print_matrix("CurEigVector", L, 1, CurEigVector, 1);

	Ap.assign(M, A);
	
	// Инициализация потоков вывода
	const int PRECISION = 10;
	vector<ofstream> fouts(L + 1);
	for (int i = 0; i < L; ++i) {
		string filename = "results/P" + to_string(i + 1) + "(t)_C++.txt";
		fouts[i].open(filename);
		fouts[i].precision(PRECISION);
	}
	fouts[L].open("results/t_C++.txt");
	fouts[L].precision(PRECISION);

	// ltmp2 == H01 + H02 + HInteration
	matadd(H01, H02, ltmp1, L);
	matadd(ltmp1, HInteration, ltmp2, L);

	for (int step = 0; step < StepsNumber; step++) {
		int dstep; // ??? Номер периода ...
		TCurrent = step * dt;
		dstep = round(TCurrent / d2);

		if (dstep < M) { // Воздействуем только M первых периодов
			At = Ap[dstep];
			Ka = 0;
			if ((TCurrent - dstep * d2 >= 0) && (TCurrent - dstep * d2 <= tau)) {
				Ka = At;
			}
			if (dstep % 2 != 0) {
				Ka = -Ka;
			}
		}
		else {
			Ka = 0;
		}

		//  Аппроксимация Паде
		vsMul(V2, Ka, ltmp3);
		matadd(ltmp2, ltmp3, Hr, L);

		vsMul(Hr, complex<TYPE>( - Ic * ddt / hh), ltmp1);
		matadd(lIdentity, ltmp1, Rr, L);
		vsMul(Hr, complex<TYPE>( + Ic * ddt / hh), ltmp1);
		matadd(lIdentity, ltmp1, Rl, L);

		linalg::inverse(Rl, L); // Вычисление обратной матрицы для Rl
	     
		dU = linalg::matmul(Rr, Rl, L, L, L, L, L, L);

		mvMul(dU, CurEigVector, UpdatedVector);
		std::swap(CurEigVector, UpdatedVector);
		// linalg::print_matrix("CurEigVector", L, 1, CurEigVector, 1);

		fouts[L] << TCurrent / 1e-9 << ' ' << Ka << '\n';
		if (Ka != 0) {
			for (int i = 0; i < L; ++i) {
				complex<TYPE> dot_product = 0;
				for (int j = 0; j < L; ++j) {
					dot_product += conj(EigVectorsR[j * L + IndexEigValuesAndVectors[i]]) * CurEigVector[j];
				}
				TYPE res = norm(dot_product);
				fouts[i] << TCurrent / 1e-9 << ' ' << res << '\n';
			}
		}
	}
	time = omp_get_wtime() - time;
	cout << "TIME ELAPSED = " << time << '\n';
	return 0;
}



