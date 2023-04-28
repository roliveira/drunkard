#ifndef UTILS_ROOT_H_
#define UTILS_ROOT_H_


#include <float.h>


// template<typename T>
// double FindLevel(T& df, double x0, double z, double eps=DBL_EPSILON) {
//     double x_old = x0;
//     double x = x0 - (IntegrateFromZeroToX(df, x0) - z) / df(x0);
//
//     while (std::abs(x - x_old) > eps) {
//         x_old = x;
//         x = x - (IntegrateFromZeroToX(df, x) - z) / df(x);
//     }
//
//     return x;
// }


template <typename F, typename T>
inline T FindInfP(F &f, T x0, int MAXITER=100) {
    int it = 0;
    T xinf = x0;
    T eps  = std::numeric_limits<T>::min();

    while (f(xinf) > eps) {
        xinf *= 10;
        ++it;
    }

    return xinf;
}

template <typename F>
inline double FindInfP(F &f, double x0, double tol, int MAXITER=100) {
    int    it   = 0;
    double xinf = x0;

    while (f(xinf) > tol) {
        xinf *= 10;
        ++it;
    }

    return xinf;
}


#endif  // UTILS_ROOT_H_
