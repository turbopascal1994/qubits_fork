#pragma once
#include <vector>
#include<vector>
#include<iostream>
#include<complex>
#include"mkl.h"
#include<iomanip>

using namespace std;

void eig(vector<complex<double>>& H0, vector<complex<double>>& EigVecLeft, vector<complex<double>>& EigVecRight, vector<complex<double>>& EigVal, int N);
lapack_int inverse(std::vector<complex<double>>& A, int n);
vector<complex<double>> getUMatrix(vector<complex<double>>& Id, vector<complex<double>>& Hr, double tstep, double h, int n);
void print_matrix(char* desc, MKL_INT m, MKL_INT n, const vector<double>& a, MKL_INT lda);
void print_matrix(char* desc, MKL_INT m, MKL_INT n, const vector<complex<double>>& a, MKL_INT lda);
vector<complex<double>> mult(vector<complex<double>>& a, vector<complex<double>>& b, int M, int N, int K, int lda, int lbd, int ldc);
vector<complex<double>> pow_matrix(vector<complex<double>> a, int step, int N);
vector<MKL_Complex16> convertToMKL(vector<complex<double>>& a);
vector<complex<double>> convertToSTD(vector<MKL_Complex16>& a);