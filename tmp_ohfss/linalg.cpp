#include "linalg.h"

void print_matrix(char* desc, MKL_INT m, MKL_INT n, const vector<double>& a, MKL_INT lda) {
    printf("\n %s\n", desc);
    for (MKL_INT i = 0; i < m; i++) {
        for (MKL_INT j = 0; j < n; j++)
            cout << fixed << setprecision(20) << ' ' << a[i * lda + j];
        cout << '\n';
    }
}

void print_matrix(char* desc, MKL_INT m, MKL_INT n, const vector<complex<double>>& a, MKL_INT lda) {
    printf("\n %s\n", desc);
    for (MKL_INT i = 0; i < m; i++) {
        for (MKL_INT j = 0; j < n; j++)
            cout << fixed << setprecision(30) << "( " << a[i * lda + j].real() << ";" << a[i * lda + j].imag() << ") ";
        cout << '\n';
    }
}

void print_eigenvalues(char* desc, MKL_INT n, double* wr, double* wi) {
    MKL_INT j;
    printf("\n %s\n", desc);
    for (j = 0; j < n; j++) {
        if (wi[j] == (double)0.0) {
            printf(" %6.2f", wr[j]);
        }
        else {
            printf(" (%6.2f,%6.2f)", wr[j], wi[j]);
        }
    }
    printf("\n");
}


void print_eigenvectors(char* desc, MKL_INT n, double* wi, double* v, MKL_INT ldv) {
    MKL_INT i, j;
    printf("\n %s\n", desc);
    for (i = 0; i < n; i++) {
        j = 0;
        while (j < n) {
            if (wi[j] == (double)0.0) {
                printf(" %6.2f", v[i + j * ldv]);
                j++;
            }
            else {
                printf(" (%6.2f,%6.2f)", v[i + j * ldv], v[i + (j + 1) * ldv]);
                printf(" (%6.2f,%6.2f)", v[i + j * ldv], -v[i + (j + 1) * ldv]);
                j += 2;
            }
        }
        printf("\n");
    }
}

vector<MKL_Complex16> convertToMKL(vector<complex<double>>& a) {
    vector<MKL_Complex16> ans(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        ans[i] = { a[i].real(), a[i].imag() };
    }
    return ans;
}

vector<complex<double>> convertToSTD(vector<MKL_Complex16>& a) {
    vector<complex<double>> ans(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        ans[i] = { a[i].real, a[i].imag };
    }
    return ans;
}

void eig(vector<complex<double>>& H0, vector<complex<double>>& EigVec, vector<complex<double>>& EigVal, int N){
    vector<MKL_Complex16> vl(N * N), vr(N * N), matrix = convertToMKL(H0), w(N * N);

    int rev = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', N, matrix.data(), N, w.data(), vl.data(), N, vr.data(), N);
    if (rev > 0) {
        cout << "Error occured after LAPACKE_zgeev\n";
        exit(1);
    }
    EigVal = convertToSTD(w);
    EigVec = convertToSTD(vl);
}   

lapack_int inverse(vector<complex<double>>& A, int n) {
    vector<lapack_int> ipiv(n + 1);
    lapack_int ret;
    auto matrix = convertToMKL(A);
    ret = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, matrix.data(), n, ipiv.data());

    if (ret != 0) return ret;

    ret = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n,  matrix.data(), n, ipiv.data());
    
    A = convertToSTD(matrix);
    return ret;
}

vector<complex<double>> mult(vector<complex<double>>& a, vector<complex<double>>& b, int M, int N, int K, int lda, int ldb, int ldc) {
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

vector<complex<double>> pow_matrix(vector<complex<double>> a, int step, int N){
    vector<complex<double>> res(N * N);
    for (size_t i = 0; i < N; ++i) res[i * N + i] = { 1, 0 };
    for(;step;step >>= 1, a = mult(a, a, N, N, N, N, N, N)){
        if (step & 1){
            res = mult(res, a, N, N, N, N, N, N);
        }
    }
    return res;
}

vector<complex<double>> getUMatrix(vector<complex<double>> & Id, vector<complex<double>>& Hr, double tstep, double h, int n){
    vector<complex<double>> up(Id.size());
    vector<complex<double>> down(Id.size());
    complex<double> I(0, 1);
    for(size_t i = 0;i < Id.size();++i){
        up[i] = Id[i] - I * Hr[i] * tstep / (h * 2);
        down[i] = Id[i] + I * Hr[i] * tstep / (h * 2);
    }

    /*print_matrix((char *)"UP", n, n, up, n);
    print_matrix((char *)"DOWN", n, n, down, n);*/

    int ret = inverse(down, n);
    if (ret != 0) exit(1);
    return mult(up, down, n, n, n, n, n, n);
}



