#pragma once
#include <limits>
using namespace std;

const double EPS = numeric_limits<double>::epsilon() * 64; // точность вычислений


bool Equal(double a, double b, double eps = EPS);
bool Not_equal(double a, double b, double eps = EPS);
bool Less(double a, double b, double eps = EPS);
bool Less_equal(double a, double b, double eps = EPS);
bool More(double a, double b, double eps = EPS);
bool More_equal(double a, double b, double eps = EPS);
