#ifndef UTILS_INVLAPLACE_HPP_
#define UTILS_INVLAPLACE_HPP_

#include <math.h>
#include <vector>


class Laplace { 
public:

    int stehfest_steps = 18;
    
    double stehfest_coeff[18] = {
        +4.960317460317460E-05, -6.095734126984128E-01, +2.745940476190476E+02, 
        -2.630695674603174E+04, +9.572572013888889E+05, -1.735869484583333E+07,
        +1.824212226472222E+08, -1.218533288309127E+09, +5.491680025283035E+09, 
        -1.736213111520684E+10, +3.945509690352738E+10, -6.526651698517500E+10, 
        +7.873006832822083E+10, -6.855644419612083E+10, +4.198434347505357E+10, 
        -1.716093471183929E+10, +4.204550039102679E+09, -4.671722265669643E+08
    };

    double ln2 = std::log(2);

    std::vector<double> StehfestCoefficients(int N);
    double Factorial(int N);
    double InverseTransform(double (*f)(double), double t) ;
};

std::vector<double> Laplace::StehfestCoefficients(int N) { 
    int N2 = static_cast<int>(N/2); 
    int NV = 2*N2; 

    std::vector<double> H = std::vector<double>(NV);

    int sign = 1; 

    if ((N2%2) != 0) {
        sign = -1;
    }

    for (int i=0; i<NV; ++i) { 
        int kmin = static_cast<int>((i+2)/2);
        int kmax = i+1; 

        if (kmax > N2) {
            kmax = N2;
        }
            
        H[i] = 0;
        sign = -sign;

        for (int k=kmin; k<=kmax; ++k) { 
            H[i] += (std::pow(k, N2)                                \
                / this->Factorial(k))*(this->Factorial(2*k)         \
                / this->Factorial(2*k-i-1)) / this->Factorial(N2-k) \
                / this->Factorial(k-1) / this->Factorial(i+1-k);
        } 

        H[i] *= sign; 
    }

    return H;
}

double Laplace::Factorial(int N) { 
    double x = 1; 

    if (N>1) { 
        for (int i=2; i<=N; ++i) {
            x = i*x; 
        }
    } 

    return x; 
} 

double Laplace::InverseTransform(double (*f)(double), double t) { 
    double ln2t = this->ln2 / t; 
    double x = 0; 
    double y = 0; 

    for (int i=0; i<stehfest_steps; ++i) { 
        x += ln2t;
        y += stehfest_coeff[i]*f(x);
    } 

    return ln2t*y; 
} 


#endif  // UTILS_INVLAPLACE_HPP_