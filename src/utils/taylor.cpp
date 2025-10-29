#include <cmath>
#include <stdexcept>

// assume x^2 + 4 for now
double f(double x) {
    return x * x + 4;
}

double secant_slope(double (*fn)(double), double x, double h) {
    /*
        MVT: Assumes f(x) differentiable for some z in [x, x + h]
        For a really small h, it sort of approximates tangent slope

        f'(z) = (f(x + h) - f(x))/h
    */
    if (h == 0.0) throw std::invalid_argument("h must be nonzero");
    return (fn(x + h) - f(x)) / h;
}

double linear_approx(double (*fn)(double), double x, double h, double z) {
    /*
        f(x + h) = f(x) + h * f'(z)
    */
    return fn(x) + h * secant_slope(fn, x, h);
}