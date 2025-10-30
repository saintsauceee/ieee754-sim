#pragma once
#include <stdexcept>

namespace num_approx {
    using Fn = double (*)(double);
    
    double secant_slope(Fn fn, double x, double h, bool central_diff=true);
    double linear_approx(Fn fn, double x, double h);
    double nth_derivative(double (*fn)(double), double x, double h, int n);
    double taylor_approx(double (*fn)(double), double x, double h, int n);
}