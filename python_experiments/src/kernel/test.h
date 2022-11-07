#include<iostream>
#include<vector>
#include<complex>
#include "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include/mkl.h"

class TestCPP {
private:
    double h, w01, w12;
    std::vector<int> a;
    std::vector<std::complex<double>> H0, EigVec, EigVal;

    std::vector<MKL_Complex16> convertToMKL(std::vector<std::complex<double>>& a) {
        std::vector<MKL_Complex16> ans(a.size());
        for (size_t i = 0; i < a.size(); i++) {
            ans[i] = { a[i].real(), a[i].imag() };
        }
        return ans;
    }

    // vector<complex<double>> convertToSTD(vector<MKL_Complex16>& a) {
    //     std::vector<complex<double>> ans(a.size());
    //     for (size_t i = 0; i < a.size(); i++) {
    //         ans[i] = { a[i].real, a[i].imag };
    //     }
    //     return ans;
    // }
    
public:
    TestCPP(double h, double w01, double w12){
        H0 = { {0, 0}, {0, 0}, {0, 0}, {0, 0}, {h * w01, 0}, {0, 0}, {0, 0}, {0, 0}, {h * (w01 + w12), 0} };
    }
    void eig(int N){
        std::cout << N << std::endl;
        std::vector<MKL_Complex16> vl(N * N), vr(N * N), matrix = convertToMKL(H0), w(N * N);

        std::cout << LAPACK_ROW_MAJOR << std::endl;
        std::cout << matrix.size() << ' ' << w.size() << ' ' << vl.size() << ' ' << vr.size() << std::endl;
        int rev = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', N, matrix.data(), N, w.data(), vl.data(), N, vr.data(), N);
        // std::vector<std::complex<double>> res(3 * 3);
        // MKL_Complex16 one = { 1, 0 };
        // MKL_Complex16 zero = { 0, 0 };
        // auto cp = H0;
        // std::cout << "JUST BEFORE ZGEMM" << std::endl;
        // cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &one, H0.data(), 1, cp.data(), 1, &zero, res.data(), 1);
        // if (rev > 0) {
        //     cout << "Error occured after LAPACKE_zgeev\n";
        //     exit(1);
        // }
        // EigVal = convertToSTD(w);
        // EigVec = convertToSTD(vl);
    }
    std::vector<std::complex<double>> get(){
        return H0;
    }
};
