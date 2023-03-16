#include "linalg.h"

using namespace std;

vector<MKL_Complex8> linalg::convertToMKL(vector<complex<float>>& a) {
    vector<MKL_Complex8> ans(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        ans[i] = { a[i].real(), a[i].imag() };
    }
    return ans;
}
vector<MKL_Complex16> linalg::convertToMKL(vector<complex<double>>& a) {
    vector<MKL_Complex16> ans(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        ans[i] = { a[i].real(), a[i].imag() };
    }
    return ans;
}

vector<complex<float>> linalg::convertToSTD(vector<MKL_Complex8>& a) {
    vector<complex<float>> ans(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        ans[i] = { a[i].real, a[i].imag };
    }
    return ans;
}
vector<complex<double>> linalg::convertToSTD(vector<MKL_Complex16>& a) {
    vector<complex<double>> ans(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        ans[i] = { a[i].real, a[i].imag };
    }
    return ans;
}

void linalg::eig(vector<complex<double>>& H0, vector<complex<double>>& EigVecLeft, vector<complex<double>>& EigVecRight, vector<complex<double>>& EigVal, int N){
    vector<MKL_Complex16> vl(N * N), vr(N * N), matrix = convertToMKL(H0), w(N * N);

    int rev = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', N, matrix.data(), N, w.data(), vl.data(), N, vr.data(), N);
    EigVal = convertToSTD(w);
    EigVecLeft = convertToSTD(vl);
    EigVecRight = convertToSTD(vr);
}   
void linalg::eig(vector<complex<float>>& H0, vector<complex<float>>& EigVecLeft, vector<complex<float>>& EigVecRight, vector<complex<float>>& EigVal, int N) {
    vector<MKL_Complex8> vl(N * N), vr(N * N), matrix = convertToMKL(H0), w(N * N);

    int rev = LAPACKE_cgeev(LAPACK_ROW_MAJOR, 'V', 'V', N, matrix.data(), N, w.data(), vl.data(), N, vr.data(), N);
    EigVal = convertToSTD(w);
    EigVecLeft = convertToSTD(vl);
    EigVecRight = convertToSTD(vr);
}

void linalg::inverse(vector<complex<double>>& A, int n) {
    vector<lapack_int> ipiv(n + 1);
    auto matrix = convertToMKL(A);
    LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, matrix.data(), n, ipiv.data());
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n,  matrix.data(), n, ipiv.data());
    A = convertToSTD(matrix);
}
void linalg::inverse(vector<complex<float>>& A, int n) {
    vector<lapack_int> ipiv(n + 1);
    auto matrix = convertToMKL(A);
    LAPACKE_cgetrf(LAPACK_ROW_MAJOR, n, n, matrix.data(), n, ipiv.data());
    LAPACKE_cgetri(LAPACK_ROW_MAJOR, n, matrix.data(), n, ipiv.data());
    A = convertToSTD(matrix);
}

vector<complex<float>> linalg::matmul(vector<complex<float>>& a, vector<complex<float>>& b, int M, int N, int K, int lda, int ldb, int ldc) {
    vector<complex<float>> res(M * N);
    MKL_Complex8 one = { 1, 0 };
    MKL_Complex8 zero = { 0, 0 };
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, &one, a.data(), lda, b.data(), ldb, &zero, res.data(), ldc);
    /*complex<float> ans;
    for(size_t i = 0;i < M;++i) {
        for (size_t j = 0; j < N; ++j) {
            ans = { 0, 0 };
            for (size_t k = 0; k < K; ++k) {
                ans += a[i * K + k] * b[k * N + j];
            }
            res[i * N + j] = ans;
        }
    }*/
    return res;
}
vector<complex<double>> linalg::matmul(vector<complex<double>>& a, vector<complex<double>>& b, int M, int N, int K, int lda, int ldb, int ldc) {
    vector<complex<double>> res(M * N);
    MKL_Complex16 one = { 1, 0 };
    MKL_Complex16 zero = { 0, 0 };
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, &one, a.data(), lda, b.data(), ldb, &zero, res.data(), ldc);
    /*complex<double> ans;
    for(size_t i = 0;i < M;++i) {
        for (size_t j = 0; j < N; ++j) {
            ans = { 0, 0 };
            for (size_t k = 0; k < K; ++k) {
                ans += a[i * K + k] * b[k * N + j];
            }
            res[i * N + j] = ans;
        }
    }*/
    return res;
}

vector<complex<double>> linalg::solve_equations(vector<complex<double>>& A, vector<complex<double>>& b, int n) {
    vector<lapack_int> ipiv(n + 1);
    auto matrix = convertToMKL(A);
    auto bb = convertToMKL(b);
    int ret = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, matrix.data(), n, ipiv.data());
    if (ret != 0) {
        cout << "error in zgetrf";
        exit(1);
    }
    ret = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', n, 1, matrix.data(), n, ipiv.data(), bb.data(), n);
    if (ret != 0) {
        cout << "error in zgetrs";
        exit(1);
    }
    return convertToSTD(bb);
}