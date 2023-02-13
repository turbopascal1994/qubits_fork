#include <vector>
#include<vector>
#include<iostream>
#include<complex>
#include"mkl.h"
#include<iomanip>

namespace linalg {
    using namespace std;

    template<typename T>
    void print_matrix(const string& s, int m, int n, const vector<T>& a, int lda) {
        cout << "\n " << s << '\n';
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                cout << fixed << setprecision(20) << ' ' << a[i * lda + j];
            cout << '\n';
        }
    }

    template<typename T>
    void print_matrix(const string& s, int m, int n, const vector<complex<T>>& a, int lda) {
        cout << "\n " << s << '\n';
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                cout << "( " << a[i * lda + j].real() << ";" << a[i * lda + j].imag() << ") ";
            cout << '\n';
        }
    }

    vector<MKL_Complex16> convertToMKL(vector<complex<double>>& a);
    vector<MKL_Complex8> convertToMKL(vector<complex<float>>& a);

    vector<complex<double>> convertToSTD(vector<MKL_Complex16>& a);
    vector<complex<float>> convertToSTD(vector<MKL_Complex8>& a);

    void eig(vector<complex<double>>& H0, vector<complex<double>>& EigVecLeft, vector<complex<double>>& EigVecRight, vector<complex<double>>& EigVal, int N);
    void eig(vector<complex<float>>& H0, vector<complex<float>>& EigVecLeft, vector<complex<float>>& EigVecRight, vector<complex<float>>& EigVal, int N);

    void inverse(vector<complex<double>>& A, int n);
    void inverse(vector<complex<float>>& A, int n);

    vector<complex<double>> solve_equations(vector<complex<double>>& A, vector<complex<double>>& b, int n);

    vector<complex<float>> matmul(vector<complex<float>>& a, vector<complex<float>>& b, int M, int N, int K, int lda, int ldb, int ldc);
    vector<complex<double>> matmul(vector<complex<double>>& a, vector<complex<double>>& b, int M, int N, int K, int lda, int ldb, int ldc);

    template<typename T>
    vector<complex<T>> matpow(vector<complex<T>> a, int step, int N) {
        vector<complex<T>> res(N * N);
        for (size_t i = 0; i < N; ++i) res[i * N + i] = { 1, 0 };
        for (; step; step >>= 1, a = matmul(a, a, N, N, N, N, N, N)) {
            if (step & 1) {
                res = matmul(res, a, N, N, N, N, N, N);
            }
        }
        return res;
    }
    
    template<typename T>
    vector<complex<T>> getUMatrix(vector<complex<T>>& Id, vector<complex<T>>& Hr, T tstep, T h, int n) {
        vector<complex<T>> up(Id.size());
        vector<complex<T>> down(Id.size());
        complex<T> I(0, 1);
        for (size_t i = 0; i < Id.size(); ++i) {
            up[i] = Id[i] - I * Hr[i] * tstep / (h * 2);
            down[i] = Id[i] + I * Hr[i] * tstep / (h * 2);
        }
        inverse(down, n);
        return matmul(up, down, n, n, n, n, n, n);
    }    
};
