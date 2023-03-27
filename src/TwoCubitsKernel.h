#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include "linalg.h"

using namespace std;

template<typename T>
struct TwoCubitsKernel {
	const T PI = acosl(-1.0);
	const T hh = 1.054 * 1e-34; // приведённая постоянная Планка
	const T F0 = 2.06 * 1e-15; // квант потока
	const T C1 = 1.0 * 1e-12; // емкость кубита
	const complex<T> Ic = { 0.0, 1.0 };
	const int Nm = 3, L = 9;

	T Theta; // Измеренный угол поворота
	T w0q1, mu1;
	T w0q2, mu2;
	T g; // параметр взаимодействия между кубитами
	T Tmax, dt, ddt;
	int StepsNumber; // Число шагов интегрирования
	int M; // Максимальное число периодов внешнего управляющего поля
	T Acon, Atar;
	T d, d2; // Период внешнего управляющего поля
	T wt; // Частота внешнего управляющего поля
	int id2, itau;
	T tau; // Длительность импульса

	vector<complex<T>> H1, H2, H1H2, NN, Identity, lIdentity, nIdentity;
	vector<complex<T>> H0q1; // Независимые гамильтонианы первого и второго кубитов
	vector<complex<T>> H0q2;

	vector<complex<T>> H12; // служебный
	vector<complex<T>> HInteration; // Гамильтониан взаимодействия
	vector<complex<T>> H01; // Гамильтониан первого кубита в контексте общей системы
	vector<complex<T>> H02; // Гамильтониан второго кубита в контексте общей системы
	vector<complex<T>> H0; // Cтационарный гамильтониан двухкубитной системы
	vector<complex<T>> EigVectorsL; // Левые собственные вектора гамильтониана двухкубитной системы
	vector<complex<T>> EigVectorsR; // Правые собственные вектора гамильтониана двухкубитной системы
	vector<complex<T>> EigValues; // Собственные числа гамильтониана двухкубитной системы
	vector<complex<T>> ltmp1, ltmp2, ltmp3;

	vector<int> IndexEigValuesAndVectors; // Перестановка номеров собственных чисел и векторов для упорядочивания их по возрастанию действительной части собственных чисел
	vector<complex<T>> V1; // Сумма операторов рождения и уничтожения для первого кубита в контексте общей системы
	vector<complex<T>> V2; // Сумма операторов рождения и уничтожения для второго кубита в контексте общей системы
	vector<T> Ta, Ca;
	vector<complex<T>> Hr; // Гамильтониан двухкубитной системы на каждом шаге по времени
	vector<complex<T>> Rr;
	vector<complex<T>> Rl;
	vector<complex<T>> dU; // Оператор системы на шаге dt
	vector<T> Ap, Aptar; // Амплитуды ...
public:
	TwoCubitsKernel(T Theta, T w0q1, T mu1, T w0q2, T mu2, T g, T wt, T tau, T Tmax, T dt, int M): 
		Theta(Theta), w0q1(w0q1), mu1(mu1), w0q2(w0q2), mu2(mu2), g(g),
		wt(wt), d(2 * PI / wt), tau(tau),
		H1(Nm * Nm), H2(Nm * Nm), H1H2(Nm * Nm), NN(Nm * Nm),
		Identity(Nm * Nm), lIdentity(L * L), nIdentity(Nm * Nm),
		H0q1(Nm * Nm), H0q2(Nm * Nm), H12(L * L), HInteration(L * L),
		H01(L * L), H02(L * L), H0(L * L), 
		EigVectorsL(L * L), EigVectorsR(L * L), EigValues(L * L),
		ltmp1(L * L), ltmp2(L * L), ltmp3(L * L),
		IndexEigValuesAndVectors(L), V1(L * L), V2(L * L),
		Tmax(Tmax), dt(dt), ddt(dt / 2), StepsNumber(round(Tmax / dt)),
		Hr(L * L), Rr(L * L), Rl(L * L), dU(L * L),
		M(M)
	{
		T V = F0 / tau;
		T Cc = Theta / (sqrtl(2 * wt * F0 * F0) / (sqrtl(hh) * sqrtl(C1)));
		Acon = 2 * Cc * V * sqrtl((hh * wt) / (C1 * 2));
		Atar = 2 * Acon;
		d2 = d / 2.0;
		id2 = floor(d2 / dt);
		itau = floor(tau / dt);
		prepareMatrices();
	}
	vector<complex<T>> getEigVector(int index) {
		vector<complex<T>> vec(L);
		for (int i = 0; i < L; ++i) {
			vec[i] = EigVectorsR[i * L + IndexEigValuesAndVectors[index]];
		}
		return vec;
	}
	vector<complex<T>> getEigVectors() {
		return EigVectorsR;
	}
	vector<complex<T>> calculateDynamic(
		vector<complex<T>> initital_state,
		T A1_stride, T A2_stride, int phi,
		bool show_progress,
		ostream& ca_out = cout,
		ostream& ta_out = cout
	) {
		prepareAmps(A1_stride, A2_stride, phi);

		// ltmp2 == H01 + H02 + HInteration
		matadd(H01, H02, ltmp1);
		matadd(ltmp1, HInteration, ltmp2);
		vector<complex<T>> updated_state(L);
		for (int step = 0; step < StepsNumber; step++) {
			// Hr = H01 + H02 + Hint + V2*Ta(jt) + V1*Ca(jt)  
			vsMul(V2, Ta[step], ltmp3);
			matadd(ltmp2, ltmp3, ltmp1);
			vsMul(V1, Ca[step], ltmp3);
			matadd(ltmp1, ltmp3, Hr);

			vsMul(Hr, complex<T>(-Ic * ddt / hh), ltmp1);
			matadd(lIdentity, ltmp1, Rr);
			vsMul(Hr, complex<T>(+Ic * ddt / hh), ltmp1);
			matadd(lIdentity, ltmp1, Rl);

			linalg::inverse(Rl, L);

			dU = linalg::matmul(Rr, Rl, L, L, L, L, L, L);

			mvMul(dU, initital_state, updated_state);
			std::swap(initital_state, updated_state);
			double norma = 0;
			for (int i = 0; i < L; ++i) {
				norma += norm(initital_state[i]);
			}
			norma = sqrt(norma);
			for (int i = 0; i < L; ++i) {
				initital_state[i] /= norma;
			}
			if (show_progress) {
				ca_out << step * dt / 1e-9 << ' ' << Ca[step] << '\n';
				ta_out << step * dt / 1e-9 << ' ' << Ta[step] << '\n';
			}
		}
		return initital_state;
	}
	vector<T> getPopulations(const vector<complex<T>>& final_state) {
		vector<T> pops(L);
		for (int i = 0; i < L; ++i) {
			complex<TYPE> dot_product = 0;
			for (int j = 0; j < L; ++j) {
				dot_product += conj(EigVectorsR[j * L + IndexEigValuesAndVectors[i]]) * final_state[j];
			}
			pops[i] = norm(dot_product);
		}
		return pops;
	}
private:
	void fillIdentity(vector<complex<T>>& A, int dim) {
		for (size_t i = 0; i < dim; ++i) {
			A[i * dim + i] = 1;
		}
	}
	void vsMul(const vector<complex<T>>& A, complex<T> b, vector<complex<T>>& Res) {
		for (size_t i = 0; i < A.size(); ++i) {
			Res[i] = A[i] * b;
		}
	}
	void matadd(const vector<complex<T>>& A, const vector<complex<T>>& B, vector<complex<T>>& Res) {
		for (size_t i = 0; i < A.size(); ++i) {
			Res[i] = A[i] + B[i];
		}
	}
	void mvMul(const vector<complex<T>>& A, const vector<complex<T>>& B, vector<complex<T>>& Res) {
		for (size_t i = 0; i < Res.size(); ++i) {
			Res[i] = 0;
			for (size_t j = 0; j < B.size(); ++j) {
				Res[i] += A[i * B.size() + j] * B[j];
			}
		}
	}
	void kMul(const vector<complex<T>>& A, const vector<complex<T>>& B, vector<complex<T>>& Res, int dim) {
		for (size_t i = 0; i < dim * dim; ++i) {
			for (size_t j = 0; j < dim * dim; ++j) {
				// Res[i][j] = A[i / Size][j / Size] * B[i % Size][j % Size];
				Res[i * dim * dim + j] = A[(i / dim) * dim + j / dim] * B[(i % dim) * dim + j % dim];
			}
		}
	}
	void prepareMatrices() {
		complex<T> zero = { 0, 0 };
		fill(H1.begin(), H1.end(), zero);
		fill(H2.begin(), H2.end(), zero);
		fill(H1H2.begin(), H1H2.end(), zero);
		fill(NN.begin(), NN.end(), zero);

		H1[1] = H2[Nm] = 1;
		for (int i = 1; i < Nm - 1; i++) {
			H1[i * Nm + i + 1] = { sqrt(i + 1), 0 };
			H2[(i + 1) * Nm + i] = { sqrt(i + 1), 0 };
		}

		NN = linalg::matmul(H2, H1, Nm, Nm, Nm, Nm, Nm, Nm);
		matadd(H1, H2, H1H2);

		fillIdentity(Identity, Nm);
		fillIdentity(lIdentity, L);
		vsMul(Identity, -1.0, nIdentity);

		// гамильтонианы отдельных кубитов размерностью NxN
		vector<complex<T>> tmp1(Nm * Nm), tmp2(Nm * Nm), tmp3(Nm * Nm), tmp4(Nm * Nm);
		vsMul(NN, complex<T>(hh * w0q1), tmp1);
		matadd(NN, nIdentity, tmp2);
		tmp3 = linalg::matmul(NN, tmp2, Nm, Nm, Nm, Nm, Nm, Nm);
		vsMul(tmp3, complex<T>(0.5 * mu1 * hh), tmp4);
		matadd(tmp1, tmp4, H0q1);

		vsMul(NN, complex<T>(hh * w0q2), tmp1);
		matadd(NN, nIdentity, tmp2);
		tmp3 = linalg::matmul(NN, tmp2, Nm, Nm, Nm, Nm, Nm, Nm);
		vsMul(tmp3, complex<T>(0.5 * mu2 * hh), tmp4);
		matadd(tmp1, tmp4, H0q2);

		// Гамильтониан взаимодействия
		kMul(H1H2, H1H2, H12, Nm);
		vsMul(H12, complex<T>(hh * g), HInteration);

		// Гамильтонианы кубитов в контексте общей системы
		kMul(H0q1, Identity, H01, Nm);
		kMul(Identity, H0q2, H02, Nm);

		// Cтационарный гамильтониан двухкубитной системы
		matadd(H01, H02, ltmp1);
		matadd(ltmp1, HInteration, H0);

		// Вычисление собственных чисел и векторов стационарного гамильтониана двухкубитной системы
		linalg::eig(H0, EigVectorsL, EigVectorsR, EigValues, L);

		// Сортировка собственных чисел по возрастанию действительной части
		iota(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), 0);
		sort(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), [&](int el1, int el2) {
			return EigValues[el1].real() < EigValues[el2].real();
		});

		kMul(H1H2, Identity, V1, Nm);
		kMul(Identity, H1H2, V2, Nm);
	}
	void prepareAmps(T A1_stride, T A2_stride, int phi) {
		Ap.assign(M, Acon);
		Aptar.assign(M, Atar);
		Ta.assign(StepsNumber, 0);
		Ca.assign(StepsNumber, 0);
		
		auto CalculateIndex = [&](double At, int j, int jt) {
			if (jt >= j * id2 && jt <= itau + j * id2) {
				if (j % 2 == 0) {
					return At;
				}
				return -At;
			}
			return 0.0;
		};
		for (int step = 0; step < StepsNumber; ++step) {
			T TCurrent = dt * step;
			int dstep = floor(TCurrent / d2);
			if (dstep < M) {
				Ta[step] = CalculateIndex(Ap[dstep] + A1_stride, dstep, step);
				if (step < phi) {
					Ca[step] = 0;
				}
				else if (step + phi < StepsNumber) {
					Ca[step + phi] = CalculateIndex(Aptar[dstep] + A2_stride, dstep, step);
				}
			}
			else {
				Ta[step] = Ca[step] = 0;
			}
		}
	}
};
