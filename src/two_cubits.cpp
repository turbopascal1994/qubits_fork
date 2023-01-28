#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <math.h>
#include <complex>
#include <mkl.h>
#include "linalg.h"

using namespace std;

const double PI = acos(-1.0);

//----------------- Constants ---------------------------- -
//real, parameter::pi = 3.1415926535897932384626433832795
//real, parameter::hh = 1.054 * 1E-34
//real, parameter::F0 = 2.06 * 1E-15
//real, parameter::C1 = 1. * 1E-12
double hh = 1.054 * 1E-34; // приведённая постоянная Планка
double F0 = 2.06 * 1E-15; // квант потока
double C1 = 1.0 * 1E-12; // емкость кубита
//
//
//!integer, parameter::L = 5.
//complex, parameter::Ic = (0., 1.)
complex<double> Ic(0.0, 1.0);
//
//integer, parameter::Nm = 3.
const int Nm = 3; // Число учитываемых уровней одного кубита
//integer, parameter::L = Nm * Nm
const int L = Nm * Nm; // Размерность системы в целом
//real, parameter::w0q1 = 5.12 * (2 * pi) * 1E9 !5.3463
double w0q1 = 5.12 * (2 * PI) * 1e9; // Собственная частота первого кубита
//real, parameter::w0q2 = 5.350 * (2 * pi) * 1E9  !5.1167
double w0q2 = 5.350 * (2 * PI) * 1e9; // Собственная частота второго кубита
//
//real, parameter::wr = 7 * (2 * pi) * 1E9
//real, parameter::g1 = 0.07 * (2 * pi) * 1E9
//real, parameter::g2 = 0.07 * (2 * pi) * 1E9
//real, parameter::g = 0.02 * (2 * pi) * 1E9!Abs((g1 * g2 * (w0q1 + w0q2 - 2 * wr)) / (2 * (w0q1 - wr) * (w0q2 - wr))) !!!0.01 * (2 * pi) * 1E9!0.0017 * (2 * pi) * 1E9!
double g = 0.02 * (2 * PI) * 1e9; // параметр взаимодействия между кубитами
//
//real, parameter::dw = w0q1 - w0q1  !5.1167
//real, parameter::wt = 5.1258 * (2 * pi) * 1E9 !g = 0.015 * (2 * pi) * 1E9
double wt = 5.1258 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
//
//real, parameter::mu1 = -0.353 * (2 * pi) * 1E9
//real, parameter::mu2 = -0.35 * (2 * pi) * 1E9
double mu1 = -0.353 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
double mu2 = -0.35 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита

//integer, parameter::M = 20000
const int M = 20000; // Максимальное число периодов внешнего управляющего поля
//real, parameter::d = 2 * pi / (wt)
double d = 2 * PI / (wt); // Период внешнего управляющего поля
//real, parameter::tau = 4 * 1E-12
double tau = 4 * 1e-12; // Длительность импульса
//
//real, parameter::dr = 4 * 1E-13
//real, parameter::dt = 2 * 1E-13, ddt = dt / 2.
double dt = 2 * 1e-13, ddt = dt / 2.0;
//real, parameter::tmax = 270 * 1E-9
double TMax = 270 * 1e-9; // Время интегрирования
//integer, parameter::Nsteps = tmax / dt
int StepsNumber = (int)(TMax / dt); // Число шагов интегрирования
//
//real, parameter::tetta = 0.0018
double Theta = 0.0018; // Измеренный угол поворота
//real, parameter::V = F0 / tau
double V = F0 / tau; // "Высота" импульса
//real, parameter::tcv = tau / 2.0
//
//!real, parameter::tx = 30 * 1E-9
//!real, parameter::Trise = 50 * 1E-9
//!integer, parameter::irise = Trise / dr
//!real, parameter::sigmax = Trise / 4.0
//!real, parameter::Ax = 0.05 * (2 * pi) * 1E9 !0.00221 Pi / 2 (k = 12)
//
//!real, parameter::tcr = 330 * 1E-9 !128.193 * 1E-9 !41.865 * 1E-9!128.193 * 1E-9
//!real, parameter::tc = 2 * Trise
//!integer, parameter::itcr = tcr / dr
//!!real, parameter::Trise = 15 * 1E-9
//!real, parameter::sigmacr = Trise / 3.0
//!real, parameter::Acr = 0.036 * (2 * pi) * 1E9 !0.0940 Pi / 2 (k = 12)
//!real, parameter::Acon = -0.001944 * (2 * pi) * 1E9 !0.00162 Pi / 2 (k = 12)
//
//integer, parameter::LWMAX = 1000
//
//complex, dimension(Nm, Nm) ::h1, h2, nn, h1h2
vector<complex<double>> H1(Nm * Nm);
vector<complex<double>> H2(Nm * Nm);
vector<complex<double>> H1H2(Nm * Nm);
vector<complex<double>> NN(Nm * Nm);
//complex, dimension(L, L) ::h12, Hint, H01, H02, H0, V1, V2, Hr, VL, VR, W, Rr, Rl, Edd, Ud, dU
// Гамильтонианы для общей системы
vector<complex<double>> H12(L * L); // служебный
vector<complex<double>> HInteration(L * L); // Гамильтониан взаимодействия
vector<complex<double>> H01(L * L); // Гамильтониан первого кубита в контексте общей системы
vector<complex<double>> H02(L * L); // Гамильтониан второго кубита в контексте общей системы
vector<complex<double>> H0(L * L); // Cтационарный гамильтониан двухкубитной системы
vector<complex<double>> EigVectorsL(L * L); // Левые собственные вектора гамильтониана двухкубитной системы
vector<complex<double>> EigVectorsR(L * L); // Правые собственные вектора гамильтониана двухкубитной системы
vector<complex<double>> EigValues(L); // Собственные числа гамильтониана двухкубитной системы
vector<int> IndexEigValuesAndVectors(L); // Перестановка номеров собственных чисел и векторов для упоредочивания их по возрастанию действительной части собственных чисел
vector<complex<double>> V1(L * L); // Сумма операторов рождения и уничтожения для первого кубита в контексте общей системы
vector<complex<double>> V2(L * L); // Сумма операторов рождения и уничтожения для второго кубита в контексте общей системы
vector<complex<double>> Hr(L * L); // Гамильтониан двухкубитной системы на каждом шаге по времени
vector<complex<double>> Rr(L * L);
vector<complex<double>> Rl(L * L);
vector<complex<double>> dU(L * L); // Оператор системы на шаге dt
//complex, dimension(Nm, Nm) ::Ed, H0q1, H0q2
vector<complex<double>> Identity(Nm * Nm);
vector<complex<double>> nIdentity(Nm * Nm);
vector<complex<double>> lIdentity(L * L);
vector<complex<double>> H0q1(Nm * Nm); // Независимые гамильтонианы первого и второго кубитов
vector<complex<double>> H0q2(Nm * Nm);
//real::buf, a11, a12, a13, a14, a15, a16, a17, a18, a19, d2
double d2; // ??? Половина периода ...
//real::Cc, A, Ka, At
double Cc; //
double A; //
double At; //
double Ka; //
//complex, dimension(L, L) ::evea, eveb
//complex, dimension(L) ::eign, fi
int CurEigVectorNumber;
vector<complex<double>> CurEigVector(L); // Собственный вектор гамильтониана двухкубитной системы, для которого выполняется интегрирование
vector<complex<double>> UpdatedVector(L);
//integer, dimension(L) ::ipiv
//integer, dimension(M + 1) ::Am
//!real, dimension(M* Nc + 1) ::Ap
//real, dimension(M + 1) ::Ap
vector<double> Ap(M); // ??? Амплитуды ...
//integer::i, jt, j, s, jj
//real::it, t
double TCurrent;
//integer::info, LWORK
//REAL             RWORK(2 * L)
//COMPLEX          WORK(LWMAX)
//real, DIMENSION(L) ::wor

vector<complex<double>> tmp1(Nm * Nm);
vector<complex<double>> tmp2(Nm * Nm);
vector<complex<double>> tmp3(Nm * Nm);
vector<complex<double>> tmp4(Nm * Nm);

vector<complex<double>> ltmp1(L * L);
vector<complex<double>> ltmp2(L * L);
vector<complex<double>> ltmp3(L * L);

void FillIdentity(vector<complex<double>>& A) {
  fill(A.begin(), A.end(), 0);
  int dim = sqrt(A.size());
  for (size_t i = 0; i < dim; ++i) {
    A[i * dim + i] = 1;
  }
}

void vAdd(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res) {
  for (size_t i = 0; i < A.size(); ++i) {
    Res[i] = A[i] + B[i];
  }
}

void vsMul(const vector<complex<double>>& A, complex<double> b, vector<complex<double>>& Res) {
  for (size_t i = 0; i < A.size(); ++i) {
    Res[i] = A[i] * b;
  }
}

void mmMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res) {
  int dim = sqrt(A.size());
  complex<double> tmp;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      tmp = 0;
      for (int k = 0; k < dim; k++) {
        tmp += A[i * dim + k] * B[k * dim + j];
      }
      Res[i * dim + j] = tmp;
    }
  }
}

void mvMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res) {
    for (size_t i = 0; i < Res.size(); ++i) {
        Res[i] = 0;
        for (size_t j = 0; j < B.size(); ++j) {
            Res[i] += A[i * B.size() + j] * B[j];
        }
    }
}


void kMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res) {
    int dim = sqrt(sqrt(Res.size()));
    for (size_t i = 0; i < dim * dim; ++i) {
        for (size_t j = 0; j < dim * dim; ++j) {
        //      Res[i][j] = A[i / Size][j / Size] * B[i % Size][j % Size];
            Res[i * dim * dim + j] = A[(i / dim) * dim + j / dim] * B[ (i % dim) * dim + j % dim];
        }
    }
}

int main() {
    Cc = Theta / (sqrt((2 * wt * F0 * F0)) / (sqrt(hh) * sqrt(C1)));
    A = 2 * Cc * V * sqrt((hh * wt) / (2 * C1));
    d2 = d / 2.0;

    // операторы рождения, уничтожения и их сумма, оператор числа частиц
    complex<double> zero = {0, 0};
    fill(H1.begin(), H1.end(), zero);
    fill(H2.begin(), H2.end(), zero);
    fill(H1H2.begin(), H1H2.end(), zero);
    fill(NN.begin(), NN.end(), zero);
    
    for (int i = 0; i < Nm - 1; i++) {
        H1[i * Nm + i + 1] = sqrt(i);
        H2[(i + 1) * Nm + i] = sqrt(i);
        H1H2[i * L + i + 1] = sqrt(i);
        H1H2[(i + 1) * L + i] = sqrt(i);
        NN[(i + 1) * Nm + i + 1] = i;
    }

    FillIdentity(Identity);
    vsMul(Identity, -1.0, nIdentity);

    // гамильтонианы отдельных кубитов размерностью NxN
    vsMul(NN, complex<double>(hh * w0q1), tmp1);
    vAdd(NN, nIdentity, tmp2);
    // tmp3 = mult(NN, tmp2, Nm, Nm, Nm, Nm, Nm, Nm);
    mmMul(NN, tmp2, tmp3);
    vsMul(tmp3, complex<double>(0.5 * mu1 * hh), tmp4);
    vAdd(tmp1, tmp4, H0q1);


    vsMul(NN, complex<double>(hh * w0q2), tmp1);
    vAdd(NN, nIdentity, tmp2);
    mmMul(NN, tmp2, tmp3);
    vsMul(tmp3, complex<double>(0.5 * mu2 * hh), tmp4);
    vAdd(tmp1, tmp4, H0q2);


    // // Гамильтониан взаимодействия
    kMul(H1H2, H1H2, H12);
    vsMul(H12, complex<double>(hh * g), HInteration);
    
    // // Гамильтонианы кубитов в контексте общей системы
    kMul(H0q1, Identity, H01);
    kMul(Identity, H0q2, H02);

    // // Cтационарный гамильтониан двухкубитной системы
    vAdd(H01, H02, ltmp1);
    vAdd(ltmp1, HInteration, H0);
    // // Вычисление собственных чисел и векторов стационарного гамильтониана двухкубитной системы
    eig(H0, EigVectorsL, EigVectorsR, EigValues, L);
    std::cout << "DONE\n";
    return 0;
    // LWORK = -1
    // call cgeev('Vectors', 'Vectors', L, H0, L, W, VL, L,VR, L, WORK, LWORK, RWORK, INFO)
    // LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    // call cgeev('Vectors', 'Vectors', L, H0, L, W, VL, L,VR, L, WORK, LWORK, RWORK, INFO)

    // // Сортировка собственных чисел по возрастанию действительной части
    iota(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), 0);
    
    sort(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), [&](int el1, int el2){
        return EigValues[el1].real() < EigValues[el2].real();
    });    

    kMul(H1H2, Identity, V1);
    kMul(Identity, H1H2, V2);

    CurEigVectorNumber = 3;
    for (int i = 0; i < Nm; i++) {
        CurEigVector[i] = EigVectorsR[i * L + IndexEigValuesAndVectors[CurEigVectorNumber]];
    }

    for (int i = 0; i < M; i++) {
        Ap[i] = A;
    }
    std::cout << "DONE\n";
    return 0;
    // for (int step = 0; step < StepsNumber; step++) {
    //     int dstep; // ??? Номер периода ...
    //     TCurrent = step * dt;
    //     dstep = floor(TCurrent / d2);

    //     if (dstep < M) { // Воздействуем только M первых периодов
    //         At = Ap[dstep];
    //         Ka = 0;
    //         if ((TCurrent - dstep * d2 >= 0) && (TCurrent - dstep * d2 <= tau)) {
    //             Ka = At;
    //         }
    //         if (dstep % 2 != 0) {
    //             Ka = -Ka;
    //         }
    //     }

    //     //  Аппроксимация Паде
    //     vAdd((complex<double>*) H01, (complex<double>*) H02, (complex<double>*) ltmp1, L);
    //     vAdd((complex<double>*) ltmp1, (complex<double>*) HInteration, (complex<double>*) ltmp2, L);
    //     vsMul((complex<double>*) V2, Ka, (complex<double>*) ltmp3, L);
    //     vAdd((complex<double>*) ltmp2, (complex<double>*) ltmp3, (complex<double>*) Hr, L);

    //     FillIdentity((complex<double>*) lIdentity, L);
    //     vsMul((complex<double>*) Hr, - Ic * ddt / hh, (complex<double>*) ltmp1, L);
    //     vAdd((complex<double>*) lIdentity, (complex<double>*) ltmp3, (complex<double>*) Rr, L);
    //     vsMul((complex<double>*) Hr, + Ic * ddt / hh, (complex<double>*) ltmp1, L);
    //     vAdd((complex<double>*) lIdentity, (complex<double>*) ltmp3, (complex<double>*) Rl, L);

    //     // Вычисление обратной матрицы для Rl
    //     //    call cgetrf(L, L, Rl, L, ipiv, info);
    //     //    call cgetri(L, Rl, L, ipiv, wor, L, info);

    //     mmMul((complex<double>*) Rr, (complex<double>*) Rl, (complex<double>*) dU, L);

    //     mvMul((complex<double>*) dU, (complex<double>*) CurEigVector, (complex<double>*) UpdatedVector, L);
    //     vsMul((complex<double>*) UpdatedVector, 1.0, (complex<double>*) CurEigVector, L);
    // }

    return 0;
}



