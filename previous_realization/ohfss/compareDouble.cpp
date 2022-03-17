#include"compareDouble.h"
#include <algorithm>
#include <iostream>

bool Equal(double a, double b, double eps) {
    return fabs(a - b) <= eps * (Less(a, b, eps) ? a : b);
}

bool Not_equal(double a, double b, double eps) {
    return !Equal(a, b, eps);
}

bool Less(double a, double b, double eps) {
    return a < b;
}

bool Less_equal(double a, double b, double eps) {
    return Less(a, b, eps) || Equal(a, b, eps);
}

bool More(double a, double b, double eps) {
    return a > b;
}

bool More_equal(double a, double b, double eps) {
    return !Less(a, b, eps);
}