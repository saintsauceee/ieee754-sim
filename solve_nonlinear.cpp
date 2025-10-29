#include <iostream>
#include <cmath>
#include <stdexcept>

// assume x^2 + 4 for now
double f(double x) {
    return pow(x, 2) + 4;
}

double bisection_method(double (*f)(double), double a, double b, double eps_x, double eps_f) {
    if (eps_x <= 0 || eps_f <= 0) {
        throw std::invalid_argument("eps must be > 0");
    }

    double fa = f(a), fb = f(b);
    if (fa == 0) return a;
    if (fb == 0) return b;

    if (fa * fb > 0) throw std::invalid_argument("interval must bracket a root");

    for (;;) {
        double c = (a + b) / 2;
        double fc = f(c);

        if (std::abs(fc) <= eps_f || (b - a) * 0.5 <= eps_x) return c;
        
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        }
        else {
            a = c;
            fa = fc;
        }
    }
}