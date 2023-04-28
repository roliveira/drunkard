#ifndef UTILS_GAMMA_HPP_
#define UTILS_GAMMA_HPP_

#include <cmath>
#include <complex>
#include <tgmath.h>

#include <boost/math/special_functions/gamma.hpp>


double gammau(double a, double x)  {
    /*if (std::abs(a) <= std::numeric_limits<float>::epsilon()*std::abs(0.0))
        return boost::math::tgamma(std::abs(a), x);
    else*/ if (a < 0)
        return (gammau(a+1, x)-std::pow(x, a)*std::exp(-x))/a;
    else
        return boost::math::tgamma(a, x);
}

std::complex<double> gammau(double a, std::complex<double> x) {

    if (a < 0) {
        return (gammau(a+1, x)-std::pow(x, a)*std::exp(-x))/a;
    } else {
        double eps = 1E-16;
        std::complex<double> b = 0.0;

        // series expansion for x < a+1

        if (std::abs(x) < a+1) {
            double ap = a;
            std::complex<double> del = 1.0;
            std::complex<double> sum = del;
            
            while (std::abs(std::norm(del)) > eps) {
                ++ap;
                del *= x / ap;
                sum += del;
            }

            b = sum * std::exp(-x + a*std::log(x) - std::lgamma(a+1));
            b = 1.0 - b;

        } else { // series expansion for x >= a+1
            std::complex<double> a0 = 1;
            std::complex<double> a1 = x;
            std::complex<double> b0 = 0;
            std::complex<double> b1 = a0;

            std::complex<double> fac = 1.0 / a1;
            int n = 1;
            
            std::complex<double> g = b1 * fac;
            std::complex<double> gold = b0;

            while (std::abs(std::norm(g-gold)) > eps) {
                gold = g;
                double ana = n - a;
                a0 = (a1 + a0 * ana) * fac;
                b0 = (b1 + b0 * ana) * fac;
                std::complex<double> anf = static_cast<std::complex<double>>(n)*fac;
                a1 = x * a0 + anf * a1;
                b1 = x * b0 + anf * b1;
                fac = 1.0 / a1;
                g = b1 * fac;

                ++n;
            }
            
            b = std::exp(-x + a * std::log(x) - std::lgamma(a)) * g;
        }

        return b*std::tgamma(a);
    }
}


#endif // UTILS_GAMMA_HPP_