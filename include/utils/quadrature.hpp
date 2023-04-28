#ifndef UTILS_QUADRATURE_H_
#define UTILS_QUADRATURE_H_


#include "qdt.h"


template<typename T> 
double IntegrateFromZeroToX(T& f, double x) {
    auto method = qdt::adaptive(FLT_EPSILON, qdt::rectangle());
    return method.integrate(f, DBL_EPSILON, x);
}

template<typename T>
double IntegrateFromZeroToInf(T& f) {
    auto method = qdt::adaptive(FLT_EPSILON, qdt::rectangle());
    return method.integrate(f, DBL_EPSILON, qdt::INF);
}

#endif  // UTILS_QUADRATURE_H_
