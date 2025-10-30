#include <cmath>
#include <stdexcept>

// assume x^2 + 4 for now
double f(double x) {
    return x * x + 4;
}

double secant_slope(double (*fn)(double), double x, double h) {
    /*
        MVT: Assumes f differentiable for some z in [x, x + h]
        For a really small h, it sort of approximates tangent slope

        f'(z) = (f(x + h) - f(x))/h
    */
    if (h == 0.0) throw std::invalid_argument("h must be nonzero");
    return (fn(x + h) - fn(x)) / h;
}

double linear_approx(double (*fn)(double), double x, double h, double z) {
    /*
        f(x + h) = f(x) + h * f'(z)
    */
    return fn(x) + h * secant_slope(fn, x, h);
}

double nth_derivative(double (*fn)(double), double x, double h, int n) {
    return 0;
}

double taylor_approx(double (*fn)(double), double x, double h, int n) {
    /*
        Taylor: Assumes f to have continuous derivative of order 0, 1, ..., n + 1

        f(x + h) = sum_k (f^k(x) * h^k) / k! + E_(n + 1) for k = 0 to n
        E_(n + 1) = (f^(n + 1)(z) * h^(n + 1)) / (n + 1)! for z in [x, x + h]
    */
    double ans = 0.0;
    double k_fact = 1.0;
    double h_pow = 1.0;

    for (int k = 0; k <= n; k++) {
        double dk = nth_derivative(fn, x, h, k);
        ans += dk * h_pow / k_fact;

        h_pow *= h;
        k_fact *= (k + 1);
    }
    return ans;
}