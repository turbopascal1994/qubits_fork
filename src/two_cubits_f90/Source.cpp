#include<iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <iomanip>
//#include <mkl.h>

using std::complex;
const double M_PI = acos(-1.0);

//----------------- Constants ---------------------------- -
//real, parameter::pi = 3.1415926535897932384626433832795
//real, parameter::hh = 1.054 * 1E-34
//real, parameter::F0 = 2.06 * 1E-15
//real, parameter::C1 = 1. * 1E-12
double hh = 1.054 * 1E-34; // ���������� ���������� ������
double F0 = 2.06 * 1E-15; // ����� ������
double C1 = 1.0 * 1E-12; // ������� ������
//
//
//!integer, parameter::L = 5.
//complex, parameter::Ic = (0., 1.)
complex<double> Ic(0.0, 1.0);
//
//integer, parameter::Nm = 3.
const int Nm = 3; // ����� ����������� ������� ������ ������
//integer, parameter::L = Nm * Nm
const int L = Nm * Nm; // ����������� ������� � �����
//real, parameter::w0q1 = 5.12 * (2 * pi) * 1E9 !5.3463
double w0q1 = 5.12 * (2 * M_PI) * 1e9; // ����������� ������� ������� ������
//real, parameter::w0q2 = 5.350 * (2 * pi) * 1E9  !5.1167
double w0q2 = 5.350 * (2 * M_PI) * 1e9; // ����������� ������� ������� ������
//
//real, parameter::wr = 7 * (2 * pi) * 1E9
//real, parameter::g1 = 0.07 * (2 * pi) * 1E9
//real, parameter::g2 = 0.07 * (2 * pi) * 1E9
//real, parameter::g = 0.02 * (2 * pi) * 1E9!Abs((g1 * g2 * (w0q1 + w0q2 - 2 * wr)) / (2 * (w0q1 - wr) * (w0q2 - wr))) !!!0.01 * (2 * pi) * 1E9!0.0017 * (2 * pi) * 1E9!
double g = 0.02 * (2 * M_PI) * 1e9; // �������� �������������� ����� ��������
//
//real, parameter::dw = w0q1 - w0q1  !5.1167
//real, parameter::wt = 5.1258 * (2 * pi) * 1E9 !g = 0.015 * (2 * pi) * 1E9
double wt = 5.1258 * (2 * M_PI) * 1e9; // ������� �������� ������������ ����
//
//real, parameter::mu1 = -0.353 * (2 * pi) * 1E9
//real, parameter::mu2 = -0.35 * (2 * pi) * 1E9
double mu1 = -0.353 * (2 * M_PI) * 1e9; // �������� ������������ ������� ������
double mu2 = -0.35 * (2 * M_PI) * 1e9; // �������� ������������ ������� ������

//integer, parameter::M = 20000
const int M = 20000; // ������������ ����� �������� �������� ������������ ����
//real, parameter::d = 2 * pi / (wt)
double d = 2 * M_PI / (wt); // ������ �������� ������������ ����
//real, parameter::tau = 4 * 1E-12
double tau = 4 * 1e-12; // ������������ ��������
//
//real, parameter::dr = 4 * 1E-13
//real, parameter::dt = 2 * 1E-13, ddt = dt / 2.
double dt = 2 * 1e-13, ddt = dt / 2.0;
//real, parameter::tmax = 270 * 1E-9
double TMax = 270 * 1e-9; // ����� ��������������
//integer, parameter::Nsteps = tmax / dt
int StepsNumber = (int)(TMax / dt); // ����� ����� ��������������
//
//real, parameter::tetta = 0.0018
double Theta = 0.0018; // ���������� ���� ��������
//real, parameter::V = F0 / tau
double V = F0 / tau; // "������" ��������
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
complex<double> H1[Nm][Nm], H2[Nm][Nm], H1H2[Nm][Nm], NN[Nm][Nm];
//complex, dimension(L, L) ::h12, Hint, H01, H02, H0, V1, V2, Hr, VL, VR, W, Rr, Rl, Edd, Ud, dU
// ������������� ��� ����� �������
complex<double> H12[L][L]; // ���������
complex<double> HInteration[L][L]; // ������������ ��������������
complex<double> H01[L][L]; // ������������ ������� ������ � ��������� ����� �������
complex<double> H02[L][L]; // ������������ ������� ������ � ��������� ����� �������
complex<double> H0[L][L]; // C����������� ������������ ������������ �������
complex<double> EigVectorsL[L][L]; // ����� ����������� ������� ������������� ������������ �������
complex<double> EigVectorsR[L][L]; // ������ ����������� ������� ������������� ������������ �������
complex<double> EigValues[L]; // ����������� ����� ������������� ������������ �������
int IndexEigValuesAndVectors[L]; // ������������ ������� ����������� ����� � �������� ��� �������������� �� �� ����������� �������������� ����� ����������� �����
complex<double> V1[L][L]; // ����� ���������� �������� � ����������� ��� ������� ������ � ��������� ����� �������
complex<double> V2[L][L]; // ����� ���������� �������� � ����������� ��� ������� ������ � ��������� ����� �������
complex<double> Hr[L][L]; // ������������ ������������ ������� �� ������ ���� �� �������
complex<double> Rr[L][L], Rl[L][L];
complex<double> dU[L][L]; // �������� ������� �� ���� dt
//complex, dimension(Nm, Nm) ::Ed, H0q1, H0q2
complex<double> Identity[Nm][Nm], nIdentity[Nm][Nm], lIdentity[L][L];
complex<double> H0q1[Nm][Nm], H0q2[Nm][Nm]; // ����������� ������������� ������� � ������� �������
//real::buf, a11, a12, a13, a14, a15, a16, a17, a18, a19, d2
double d2; // ??? �������� ������� ...
//real::Cc, A, Ka, At
double Cc; //
double A; //
double At; //
double Ka; //
//complex, dimension(L, L) ::evea, eveb
//complex, dimension(L) ::eign, fi
int CurEigVectorNumber;
complex<double> CurEigVector[L]; // ����������� ������ ������������� ������������ �������, ��� �������� ����������� ��������������
complex<double> UpdatedVector[L];
//integer, dimension(L) ::ipiv
//integer, dimension(M + 1) ::Am
//!real, dimension(M* Nc + 1) ::Ap
//real, dimension(M + 1) ::Ap
double Ap[M]; // ??? ��������� ...
//integer::i, jt, j, s, jj
//real::it, t
double TCurrent;
//integer::info, LWORK
//REAL             RWORK(2 * L)
//COMPLEX          WORK(LWMAX)
//real, DIMENSION(L) ::wor

complex<double> tmp1[Nm][Nm], tmp2[Nm][Nm], tmp3[Nm][Nm], tmp4[Nm][Nm];
complex<double> ltmp1[L][L], ltmp2[L][L], ltmp3[L][L];

void FillIdentity(complex<double>* A, int Size);
void vAdd(complex<double>* A, complex<double>* B, complex<double>* Res, int Size);
void vsMul(complex<double>* A, complex<double> b, complex<double>* Res, int Size);
void mmMul(complex<double>* A, complex<double>* B, complex<double>* Res, int Size);
void mvMul(complex<double>* A, complex<double>* B, complex<double>* Res, int Size);
void kMul(complex<double>* A, complex<double>* B, complex<double>* Res, int Size);
int CompareEigValues(const void* elem1, const void* elem2);

int main(int argc, char* argv[]) {
  Cc = Theta / (sqrt((2 * wt * F0 * F0)) / (sqrt(hh) * sqrt(C1)));
  A = 2 * Cc * V * sqrt((hh * wt) / (2 * C1));
  d2 = d / 2.0;

  // ��������� ��������, ����������� � �� �����, �������� ����� ������
  for (int i = 0; i < Nm; i++) {
    for (int j = 0; j < Nm; j++) {
      H1[i][j] = 0;
      H2[i][j] = 0;
      H1H2[i][j] = 0;
      NN[i][j] = 0;
    }
  }
  for (int i = 0; i < Nm - 1; i++) {
    H1[i][i + 1] = sqrt((double)i);
    H2[i + 1][i] = sqrt((double)i);
    H1H2[i][i + 1] = sqrt((double)i);
    H1H2[i + 1][i] = sqrt((double)i);
    NN[i + 1][i + 1] = (double)i;
  }

  FillIdentity((complex<double>*) Identity, Nm);
  vsMul((complex<double>*) Identity, -1.0, (complex<double>*) nIdentity, Nm);

// ������������� ��������� ������� ������������ NxN
  vsMul((complex<double>*) NN, hh * w0q1, (complex<double>*) tmp1, Nm);
  vAdd((complex<double>*) NN, (complex<double>*) nIdentity, (complex<double>*) tmp2, Nm);
  mmMul((complex<double>*) NN, (complex<double>*) tmp2, (complex<double>*) tmp3, Nm);
  vsMul((complex<double>*) tmp3, 0.5 * mu1 * hh, (complex<double>*) tmp4, Nm);
  vAdd((complex<double>*) tmp1, (complex<double>*) tmp4, (complex<double>*) H0q1, Nm);


  vsMul((complex<double>*) NN, hh * w0q2, (complex<double>*) tmp1, Nm);
  vAdd((complex<double>*) NN, (complex<double>*) nIdentity, (complex<double>*) tmp2, Nm);
  mmMul((complex<double>*) NN, (complex<double>*) tmp2, (complex<double>*) tmp3, Nm);
  vsMul((complex<double>*) tmp3, 0.5 * mu2 * hh, (complex<double>*) tmp4, Nm);
  vAdd((complex<double>*) tmp1, (complex<double>*) tmp4, (complex<double>*) H0q2, Nm);

  // ������������ ��������������
  kMul((complex<double>*) H1H2, (complex<double>*) H1H2, (complex<double>*) H12, Nm);
  vsMul((complex<double>*) H12, hh * g, (complex<double>*) HInteration, Nm);
  for(int i = 0;i < L;++i){
    for(int j = 0;j < L;++j) {
      std::cout << std::fixed << std::setprecision(5) << HInteration[i][j] << ' ';
    }
    std::cout << '\n';
  }
  return 0;

  // ������������� ������� � ��������� ����� �������
  kMul((complex<double>*) H0q1, (complex<double>*) Identity, (complex<double>*) H01, Nm);
  kMul((complex<double>*) Identity, (complex<double>*) H0q2, (complex<double>*) H02, Nm);

  // C����������� ������������ ������������ �������
  vAdd((complex<double>*) H01, (complex<double>*) H02, (complex<double>*) ltmp1, L);
  vAdd((complex<double>*) ltmp1, (complex<double>*) HInteration, (complex<double>*) H0, L);

  // ���������� ����������� ����� � �������� ������������� ������������� ������������ �������
//LWORK = -1
//call cgeev('Vectors', 'Vectors', L, H0, L, W, VL, L,VR, L, WORK, LWORK, RWORK, INFO)
//LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
//call cgeev('Vectors', 'Vectors', L, H0, L, W, VL, L,VR, L, WORK, LWORK, RWORK, INFO)
  
  // ���������� ����������� ����� �� ����������� �������������� �����
  for (int i = 0; i < L; i++) {
    IndexEigValuesAndVectors[i] = i;
  }
  std::qsort(IndexEigValuesAndVectors, L, sizeof(int), CompareEigValues);

  kMul((complex<double>*) H1H2, (complex<double>*) Identity, (complex<double>*) V1, Nm);
  kMul((complex<double>*) Identity, (complex<double>*) H1H2, (complex<double>*) V2, Nm);

  CurEigVectorNumber = 3;
  for (int i = 0; i < Nm; i++) {
    CurEigVector[i] = EigVectorsR[i][IndexEigValuesAndVectors[CurEigVectorNumber]];
  }
  
  for (int i = 0; i < M; i++) {
    Ap[i] = A;
  }

  for (int step = 0; step < StepsNumber; step++) {
    int dstep; // ??? ����� ������� ...
    TCurrent = step * dt;
    dstep = floor(TCurrent / d2);

    if (dstep < M) { // ������������ ������ M ������ ��������
      At = Ap[dstep];
      Ka = 0;
      if ((TCurrent - dstep * d2 >= 0) && (TCurrent - dstep * d2 <= tau)) {
        Ka = At;
      }
      if (dstep % 2 != 0) {
        Ka = -Ka;
      }
    }

    //  ������������� ����
    vAdd((complex<double>*) H01, (complex<double>*) H02, (complex<double>*) ltmp1, L);
    vAdd((complex<double>*) ltmp1, (complex<double>*) HInteration, (complex<double>*) ltmp2, L);
    vsMul((complex<double>*) V2, Ka, (complex<double>*) ltmp3, L);
    vAdd((complex<double>*) ltmp2, (complex<double>*) ltmp3, (complex<double>*) Hr, L);

    FillIdentity((complex<double>*) lIdentity, L);
    vsMul((complex<double>*) Hr, - Ic * ddt / hh, (complex<double>*) ltmp1, L);
    vAdd((complex<double>*) lIdentity, (complex<double>*) ltmp3, (complex<double>*) Rr, L);
    vsMul((complex<double>*) Hr, + Ic * ddt / hh, (complex<double>*) ltmp1, L);
    vAdd((complex<double>*) lIdentity, (complex<double>*) ltmp3, (complex<double>*) Rl, L);

    // ���������� �������� ������� ��� Rl
//    call cgetrf(L, L, Rl, L, ipiv, info);
//    call cgetri(L, Rl, L, ipiv, wor, L, info);

    mmMul((complex<double>*) Rr, (complex<double>*) Rl, (complex<double>*) dU, L);

    mvMul((complex<double>*) dU, (complex<double>*) CurEigVector, (complex<double>*) UpdatedVector, L);
    vsMul((complex<double>*) UpdatedVector, 1.0, (complex<double>*) CurEigVector, L);
  }

  return 0;
}


void FillIdentity(complex<double>* A, int Size) {
  for (int i = 0; i < Size * Size; i++) {
    A[i] = 0;
  }
  for (int i = 0; i < Size; i++) {
    A[i * Size + i] = 1.0;
  }
  return;
}

void vAdd(complex<double>* A, complex<double>* B, complex<double>* Res, int Size) {
  for (int i = 0; i < Size; i++) {
    Res[i] = A[i] + B[i];
  }
}

void vsMul(complex<double>* A, complex<double> b, complex<double>* Res, int Size) {
  for (int i = 0; i < Size; i++) {
    Res[i] = A[i] * b;
  }
}

void mmMul(complex<double>* A, complex<double>* B, complex<double>* Res, int Size) {
  complex<double> tmp;
  for (int i = 0; i < Size; i++) {
    for (int j = 0; j < Size; j++) {
      tmp = 0;
      for (int k = 0; k < Size; k++) {
        tmp += A[i * Size + k] * B[k * Size + j];
      }
      Res[i * Size + j] = tmp;
    }
  }
}

void mvMul(complex<double>* A, complex<double>* B, complex<double>* Res, int Size) {
  for (int i = 0; i < Size; i++) {
    Res[i] = 0;
    for (int j = 0; j < Size; j++) {
      Res[i] += A[i * Size + j] * B[j];
    }
  }
}


void kMul(complex<double>* A, complex<double>* B, complex<double>* Res, int Size) {
  for (int i = 0; i < Size*Size; i++) {
    for (int j = 0; j < Size*Size; j++) {
//      Res[i][j] = A[i / Size][j / Size] * B[i % Size][j % Size];
      Res[i * Size * Size + j] = A[(i / Size) * Size + j / Size] * B[ (i % Size) * Size + j % Size];
    }
  }
}

int CompareEigValues(const void* elem1, const void* elem2) {
  if (EigValues[*(int*)elem1].real() < EigValues[*(int*)elem2].real())
    return -1;
  if (EigValues[*(int*)elem1].real() > EigValues[*(int*)elem2].real())
    return 1;
  return 0;
}
