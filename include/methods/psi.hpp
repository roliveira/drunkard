#ifndef METHODS_PSI_HPP_
#define METHODS_PSI_HPP_

#include <cmath>
#include <vector>
#include <map>

#include "io/echo.hpp"
#include "utils/invlaplace.hpp"
#include "utils/gamma.hpp"


// ====================================== //
//                                        //
//           Psi for using beta           //
//                                        //
// ====================================== //

class Psi
{

public:

    int _maxiter = 1E+2;
    double _tol   = 1E-6;

    double _t1;
    double _t2;
    double _beta;
    double _A;

    Psi(void);
    Psi(double t1, double t2, double beta);

    double NormalizationFactor(void);
    double Evaluate(double t);
    double CumulativeDensityFunction(double t);

    double SampleTime(double z);

};


Psi::Psi(void)
{}

Psi::Psi(double t1, double t2, double beta)
: _t1(t1), _t2(t2), _beta(beta)
{
    _A = 1.0/this->NormalizationFactor();
}

double Psi::NormalizationFactor(void)
{
    if (_t1 > _t2 * 1000)
    {
        return _t2;
    }
    else
    {
        return _t1*std::exp(_t1/_t2)*std::pow(_t1/_t2,_beta)*gammau(-_beta,_t1/_t2);
    }
}

double Psi::Evaluate(double t)
{
    if (_t1 > _t2 * 1000)
    {
        return _A*std::exp(-t/_t2);
    }
    else
    {
        return _A*std::exp(-t/_t2)*std::pow(1+t/_t1, -1-_beta);
    }
}

double Psi::CumulativeDensityFunction(double t)
{
    if (_t1 > _t2 * 1000)
    {
        return 1-std::exp(-t/_t2);
    }
    else
    {
        return 1-_A*_t1*std::exp(_t1/_t2)*std::pow(_t1/_t2,_beta)*gammau(-_beta, (t+_t1)/_t2);
    }
}

double Psi::SampleTime(double z)
{
    double xa = 0.0;
    double xb = _t2*1000.0;

    int niter = 0;
    double xc, fc;

    while (niter < _maxiter) 
    {
        xc = (xa+xb)/2.0;
        fc = this->CumulativeDensityFunction(xc);

        if (std::abs(fc-z) < _tol || (std::abs(xb-xa)/2.0 < _tol)) break;
        ++niter;

        if (fc < z) xa = xc;
        else        xb = xc;
    }
    
    return xc;
}


std::vector<Psi> CreateLocalPsi(std::vector<std::map<std::string, double>> constants)
{
    std::vector<Psi> out(constants.size());

    for (int i=0; i<constants.size(); ++i) 
    {
        out[i] = Psi(constants[i]["t1"], constants[i]["t2"], constants[i]["beta"]);
    }

    return out;

}

std::ostream &operator<<(std::ostream &os, const std::shared_ptr<Psi> &psi) 
{
    if (psi->_t1 > psi->_t2 * 1000)
    {
        os << message("Psi", "A*exp(-t/t2)");
    }
    else
    {
        os << message("Psi", "A*exp(-t/t2)*pow(1+t/t1,-1-beta)");
    }

    os << message("Adv. transit time",  psi->_t1  );
    os << message("Diff. cut-off time", psi->_t2  );
    os << message("Beta factor",        psi->_beta);
    os << message("Norm. factor",       psi->_A   );

    return os;
}

// ====================================== //
//                                        //
//     Psi for exact ADE in a network     //
//                                        //
// ====================================== //

class PsiADE
{
public:
    int _maxiter = 1E+2;
    double _tol   = 1E-6;

    int _stehfest_steps=18;
    double _stehfest_coeff[18] = {
        +4.960317460317460E-05, -6.095734126984128E-01, +2.745940476190476E+02, 
        -2.630695674603174E+04, +9.572572013888889E+05, -1.735869484583333E+07,
        +1.824212226472222E+08, -1.218533288309127E+09, +5.491680025283035E+09, 
        -1.736213111520684E+10, +3.945509690352738E+10, -6.526651698517500E+10, 
        +7.873006832822083E+10, -6.855644419612083E+10, +4.198434347505357E+10, 
        -1.716093471183929E+10, +4.204550039102679E+09, -4.671722265669643E+08
    };

    double _ln2 = std::log(2);

    int _nlinks;

    std::vector<double> _v;    //!< darcy velocity
    std::vector<double> _d;    //!< diffusion coefficient
    std::vector<double> _pe;   //!< peclet number
    std::vector<double> _l;    //!< length
    std::vector<double> _A;

    PsiADE(void);
    PsiADE(std::vector<double> v, double d, std::vector<double> pe, std::vector<double> l);
    PsiADE(std::vector<double> v, std::vector<double> d, std::vector<double> pe, std::vector<double> l);

    double Sigma(double s, double v, double d);
    double Bs(double s);
    double Fs(double s, int ilink);

    void NormalizationFactor(void);

    double InverseTransform(double t, int ilink);
    double SampleTime(double z, int inlink);

};

PsiADE::PsiADE(void)
{}

PsiADE::PsiADE(std::vector<double> v, double d, std::vector<double> pe, std::vector<double> l):
_v(v), _pe(pe), _l(l)
{
    this->_nlinks = this->_v.size();
    this->_d = std::vector<double>(this->_nlinks, d);
    this->NormalizationFactor();
}

PsiADE::PsiADE(std::vector<double> v, std::vector<double> d, std::vector<double> pe, std::vector<double> l):
_v(v), _d(d), _pe(pe), _l(l)
{
    this->_nlinks = this->_v.size();
    this->NormalizationFactor();
}

void PsiADE::NormalizationFactor(void)
{
    this->_A = std::vector<double>(this->_nlinks, 1.0);

    for (int i=0; i<this->_nlinks; ++i)
    {
        double temp = 1000*std::pow(this->_l[i], 2)/(2*this->_d[i]);
        this->_A[i] = 1.0/this->InverseTransform(temp, i);
    }
}

double PsiADE::InverseTransform(double t, int ilink) 
{ 
    double ln2t = this->_ln2/t;
    double x = 0.0;
    double y = 0.0;

    for (int i=0; i<this->_stehfest_steps; ++i) 
    { 
        x += ln2t;
        y += this->_stehfest_coeff[i]*this->Fs(x, ilink);
    } 

    return ln2t*y; 
} 

double PsiADE::Sigma(double s, double v, double d)
{
    return std::sqrt(std::pow(v,2.0)+4.0*d*s)/(2.0*d);
}

double PsiADE::Bs(double s)
{
    double b = 0.0;

    for (int i=0; i<this->_nlinks; ++i)
    {
        double sig = this->Sigma(s, this->_v[i], this->_d[i]);
        b += this->_d[i]*sig/std::tanh(sig*this->_l[i]);
    }

    return 1.0/b;    
}

double PsiADE::Fs(double s, int ilink)
{
    double sig = this->Sigma(s, this->_v[ilink], this->_d[ilink]);
    double b = this->Bs(s);
    double temp = b*this->_d[ilink]*sig/(std::sinh(sig*this->_l[ilink])*s);

    if (this->_v[ilink] > 0)
    {
        return temp*std::exp(this->_pe[ilink]/2);
    }
    else
    {
        return temp*std::exp(-this->_pe[ilink]/2);
    }
}

double PsiADE::SampleTime(double z, int ilink)
{
    double xa = 0.0;
    double xb = 1000.0*std::pow(this->_l[ilink], 2)/(2*this->_d[ilink]);

    int niter = 0;
    double xc, fc;

    while (niter < this->_maxiter) 
    {
        xc = (xa+xb)/2.0;
        fc = this->_A[ilink]*this->InverseTransform(xc, ilink);

        if (std::abs(fc-z) < this->_tol || (std::abs(xb-xa)/2.0 < this->_tol)) 
        {
            break;
        }

        ++niter;

        if (fc < z) xa = xc;
        else        xb = xc;
    }
    
    return xc;
}





// ================================================ //
//                                                  //
//                   Psi prototype                  //
//                                                  //
// ================================================ //

class PsiBase {
public:

    double ln2 = 0.69314718056;  // ln(2)
    int stehfest_steps = 18;     // 18th first stehfest coefficients 
    double stehfest_coeff[18] = {
        +4.960317460317460E-05, -6.095734126984128E-01, +2.745940476190476E+02, 
        -2.630695674603174E+04, +9.572572013888889E+05, -1.735869484583333E+07,
        +1.824212226472222E+08, -1.218533288309127E+09, +5.491680025283035E+09, 
        -1.736213111520684E+10, +3.945509690352738E+10, -6.526651698517500E+10, 
        +7.873006832822083E+10, -6.855644419612083E+10, +4.198434347505357E+10, 
        -1.716093471183929E+10, +4.204550039102679E+09, -4.671722265669643E+08
    };

    double tc;      // characteristic time
    double v;       // velocity
    double d;       // generalized dispersion coefficient
    double alpha;   // dimensionless dispersivity
    double alphai;  // 1 / alpha

    PsiBase(void);
    PsiBase(double tc, double v, double d);

    virtual std::complex<double> EvaluateFunction(std::complex<double> s) = 0;
    virtual double InvertFunction(double t);

    void CalculateDispersivity(void);
    std::complex<double> MemoryFunction(std::complex<double> s);
    std::complex<double> ConcentrationDirichlet(double x, std::complex<double> s);
    double InvertConcentrationDirichlet(double x, double t);
};

PsiBase::PsiBase(void): tc(0), v(0), d(0), alpha(0) {
}

PsiBase::PsiBase(double tc, double v, double d): tc(tc), v(v), d(d) {
    this->CalculateDispersivity();
}

void PsiBase::CalculateDispersivity(void) {
    this->alpha = this->d/this->v;
    this->alphai = 1.0 / this->alpha;
}

std::complex<double> PsiBase::MemoryFunction(std::complex<double> s) {
    return this->tc * s * this->EvaluateFunction(s) / (1.0 - this->EvaluateFunction(s));
}

double PsiBase::InvertFunction(double t) {
    // double ln2t = this->ln2 / t; 
    // double s = 0; 
    // double y = 0; 

    // for (int i=0; i<this->stehfest_steps; ++i) { 
    //     s += ln2t;
    //     y += this->stehfest_coeff[i]*this->EvaluateFunction(s);
    // } 

    // return ln2t*y;

    double tolerance = 1e-8;
    int M = 18;
    double dehoog_factor = 4.0;  //DeHoog time factor
    double T;                    //Period of DeHoog Inversion formula
    double gamma;                //Integration limit parameter
    std::complex<double> h2M, R2M, z, dz, s;   //Temporary variables

    std::vector<std::complex<double>> Fctrl(2*M+1);

    std::vector<std::complex<double>> d(2*M+1);
    std::vector<std::complex<double>> A(2*M+2);
    std::vector<std::complex<double>> B(2*M+2);

    std::vector<std::vector<std::complex<double>>> e(2*M, std::vector<std::complex<double>>(M));
    std::vector<std::vector<std::complex<double>>> q(2*M, std::vector<std::complex<double>>(M));

    //Calculate period and integration limits
    T = dehoog_factor*t;    
    gamma = -0.5*std::log(tolerance)/T;

    //Calculate F(s) at evalution points gamma+IM*i*PI/T for i=0 to 2*M-1
    //This is likely the most time consuming portion of the DeHoog algorithm
    Fctrl[0] = 0.5*this->EvaluateFunction(gamma); 
    
    for (int i = 1; i <= 2*M; ++i) {
        s = std::complex<double>(gamma, i*M_PI/T);
        Fctrl[i] = this->EvaluateFunction(s.real());
    }

    //Evaluate e and q
    //eqn 20 of De Hoog et al 1982
    for (int i = 0; i < 2*M; ++i) {
        e[i][0] = 0.0;
        q[i][1] = Fctrl[i+1] / Fctrl[i];
    }

    e[2*M-1][0] = 0.0;

    //one minor correction - does not work for r<=M, as in the paper
    for (int r = 1; r <= M-1; ++r) { 
        for (int i = 2*(M-r); i >= 0; --i) {
            if ((i<2*(M-r)) && (r>1)){
                q[i][r] = q[i+1][r-1]*e[i+1][r-1]/e[i][r-1];
            }

            e[i][r] = q[i+1][r]-q[i][r]+e[i+1][r-1];
        }
    }

    //Populate d vector 
    d[0] = Fctrl[0];

    for (int m = 1; m <= M; ++m) {
        d[2*m-1] = -q[0][m];
        d[2*m] = -e[0][m];
    }

    //Evaluate A, B
    //Eqn. 21 in De Hoog et al.
    z = std::complex<double>(std::cos(M_PI*t/T), std::sin(M_PI*t/T));

    A[0] = 0.0; 
    B[0] = 1.0; //A_{-1},B_{-1} in De Hoog
    A[1] = d[0]; 
    B[1] = 1.0;

    for (int n = 2; n <= 2*M+1; n++) {
        dz = d[n-1]*z; 
        A[n] = A[n-1] + dz*A[n-2];
        B[n] = B[n-1] + dz*B[n-2];
    }

    //Eqn. 23 in De Hoog et al.
    h2M = 0.5*(1.0+z*(d[2*M-1]-d[2*M]));
    R2M = -h2M*(1.0-std::sqrt(1.0+(z*d[2*M]/h2M/h2M)));

    //Eqn. 24 in De Hoog et al.
    A[2*M+1] = A[2*M] + R2M*A[2*M-1];
    B[2*M+1] = B[2*M] + R2M*B[2*M-1];

    //Final result: A[2*M]/B[2*M]=sum [F(gamma+itheta)*exp(itheta)]
    return 1.0/T*std::exp(gamma*t)*(A[2*M+1]/B[2*M+1]).real();
}

double PsiBase::InvertConcentrationDirichlet(double x, double t) {
    // double ln2t = this->ln2 / t; 
    // double s = 0; 
    // double y = 0; 

    // for (int i=0; i<this->stehfest_steps; ++i) { 
    //     s += ln2t;
    //     y += this->stehfest_coeff[i]*this->ConcentrationDirichlet(x, s);
    // } 

    // return ln2t*y;
    double tolerance = 1e-8;
    int M = 18;
    double dehoog_factor = 4.0;  //DeHoog time factor
    double T;                    //Period of DeHoog Inversion formula
    double gamma;                //Integration limit parameter
    std::complex<double> h2M, R2M, z, dz, s;   //Temporary variables

    std::vector<std::complex<double>> Fctrl(2*M+1);

    std::vector<std::complex<double>> d(2*M+1);
    std::vector<std::complex<double>> A(2*M+2);
    std::vector<std::complex<double>> B(2*M+2);

    std::vector<std::vector<std::complex<double>>> e(2*M, std::vector<std::complex<double>>(M));
    std::vector<std::vector<std::complex<double>>> q(2*M, std::vector<std::complex<double>>(M));

    //Calculate period and integration limits
    T = dehoog_factor*t;    
    gamma = -0.5*std::log(tolerance)/T;

    //Calculate F(s) at evalution points gamma+IM*i*PI/T for i=0 to 2*M-1
    //This is likely the most time consuming portion of the DeHoog algorithm
    Fctrl[0] = 0.5*this->ConcentrationDirichlet(x, gamma); 
    
    for (int i = 1; i <= 2*M; ++i) {
        s = std::complex<double>(gamma, i*M_PI/T);
        Fctrl[i] = this->ConcentrationDirichlet(x, s);
    }

    //Evaluate e and q
    //eqn 20 of De Hoog et al 1982
    for (int i = 0; i < 2*M; ++i) {
        e[i][0] = 0.0;
        q[i][1] = Fctrl[i+1] / Fctrl[i];
    }

    e[2*M-1][0] = 0.0;

    //one minor correction - does not work for r<=M, as in the paper
    for (int r = 1; r <= M-1; ++r) { 
        for (int i = 2*(M-r); i >= 0; --i) {
            if ((i<2*(M-r)) && (r>1)){
                q[i][r] = q[i+1][r-1]*e[i+1][r-1]/e[i][r-1];
            }

            e[i][r] = q[i+1][r]-q[i][r]+e[i+1][r-1];
        }
    }

    //Populate d vector 
    d[0] = Fctrl[0];

    for (int m = 1; m <= M; ++m) {
        d[2*m-1] = -q[0][m];
        d[2*m] = -e[0][m];
    }

    //Evaluate A, B
    //Eqn. 21 in De Hoog et al.
    z = std::complex<double>(std::cos(M_PI*t/T), std::sin(M_PI*t/T));

    A[0] = 0.0; 
    B[0] = 1.0; //A_{-1},B_{-1} in De Hoog
    A[1] = d[0]; 
    B[1] = 1.0;

    for (int n = 2; n <= 2*M+1; n++) {
        dz = d[n-1]*z; 
        A[n] = A[n-1] + dz*A[n-2];
        B[n] = B[n-1] + dz*B[n-2];
    }

    //Eqn. 23 in De Hoog et al.
    h2M = 0.5*(1.0+z*(d[2*M-1]-d[2*M]));
    R2M = -h2M*(1.0-std::sqrt(1.0+(z*d[2*M]/h2M/h2M)));

    //Eqn. 24 in De Hoog et al.
    A[2*M+1] = A[2*M] + R2M*A[2*M-1];
    B[2*M+1] = B[2*M] + R2M*B[2*M-1];

    //Final result: A[2*M]/B[2*M]=sum [F(gamma+itheta)*exp(itheta)]
    return 1.0/T*std::exp(gamma*t)*(A[2*M+1]/B[2*M+1]).real();
}

std::complex<double> PsiBase::ConcentrationDirichlet(double x, std::complex<double> s) {
    /*double w = std::sqrt(1.0 + 4.0 * s * this->alpha / (this->MemoryFunction(s) * this->v));
    double z = x * (1.0 + w) / (2.0 * this->alpha);
    double k = (2.0 * w - x * w + x) / (2.0 * this->alpha);

    double term1 = (w+1) + (w-1)*std::exp((x-1)*w/this->alpha); 
    double term2;

    if (w/this->alpha > 50) {
        term2 = (w+1);
    } else {
        term2 = (w+1)*std::exp(w/this->alpha) + (w-1);
    }

    std::cout << "x     = " << x << std::endl;
    std::cout << "s     = " << s << std::endl;
    std::cout << "a     = " << this->alpha << std::endl;
    std::cout << "m     = " << this->MemoryFunction(s) << std::endl;
    std::cout << "v     = " << this->v << std::endl;
    std::cout << "w     = " << w << std::endl;
    std::cout << "z     = " << z << std::endl;
    std::cout << "k     = " << k << std::endl;
    std::cout << "term1 = " << term1 << std::endl;
    std::cout << "term2 = " << term2 << std::endl;
    std::cout << "res   = " << std::exp(k) * term1 / (s * term2) << std::endl;

    return std::exp(k) * term1 / (s * term2);*/

    std::complex<double> z = std::sqrt(1.0 + 4.0 * s * this->alpha / (this->MemoryFunction(s) * this->v)) / this->alpha;
    std::complex<double> term1 = 2.0 * z * std::exp((this->alphai + z) / 2.0);
    std::complex<double> term2 = std::exp(z) * (z + this->alphai + 2.0 * s / (this->MemoryFunction(s) * this->v));
    std::complex<double> term3 = z - this->alphai - 2.0 * s / (this->MemoryFunction(s) * this->v);

    if (this->alphai > 10) {
        return (1.0 / s) * std::exp(this->alphai / 2.0 - std::sqrt(s * this->alphai / (this->MemoryFunction(s) * this->v) + this->alphai * this->alphai / 4.0));
    } else {
        return (1.0 / s) * term1 / (term2 + term3);
    }

}

class PsiTPL: public PsiBase {
public:
    double t1;     // advective time
    double t2;     // diffusive cut-off time
    double beta;   // heterogeneity parameter
    double tau2;   // t2/t1
    double tau2i;  // t1/t2
    double c;      // normalisation factor

    PsiTPL(void);
    PsiTPL(double t1, double t2, double beta);
    PsiTPL(double v, double d, double t1, double t2, double beta);

    std::complex<double> EvaluateFunction(std::complex<double> s);
    double InvertFunction(double t);

    void CalculateTau2(void);
    void CalculateNormFactor(void);
};

PsiTPL::PsiTPL(void): PsiBase(), t1(0), t2(0), beta(0), tau2(0), c(1) {
}

PsiTPL::PsiTPL(double t1, double t2, double beta): PsiBase(t1, 0, 0), t1(t1), t2(t2), beta(beta) {
    this->CalculateTau2();
    this->CalculateNormFactor();
}

PsiTPL::PsiTPL(double v, double d, double t1, double t2, double beta): PsiBase(t1, v, d), t1(t1), t2(t2), beta(beta) {
    this->CalculateTau2();
    this->CalculateNormFactor();
}

void PsiTPL::CalculateTau2(void) {
    this->tau2 = this->t2/this->t1;
    this->tau2i = this->t1/this->t2;
}

void PsiTPL::CalculateNormFactor(void) {
    this->c = this->t1 * std::pow(this->t1 / this->t2, this->beta) *
              std::exp(this->t1 / this->t2) * gammau(-this->beta, this->t1 / this->t2);
    this->c = 1.0 / this->c;
}

std::complex<double> PsiTPL::EvaluateFunction(std::complex<double> s) {
    return std::pow(1.0 + s * this->tau2 * this->t1, this->beta) * std::exp(s * this->t1) *
           gammau(-this->beta, this->tau2i + s * this->t1) / gammau(-this->beta, this->tau2i);
}

double PsiTPL::InvertFunction(double t) {
    return this->c * std::exp(-t / this->t2) * std::pow(1 + t / this->t1, -1 - this->beta);
}

class PsiEXP: public PsiBase {
public:
    double a;

    PsiEXP(void);
    PsiEXP(double a);

    double EvaluateFunction(double s);
};

PsiEXP::PsiEXP(void) {
}

PsiEXP::PsiEXP(double a): a(a) {
}

double PsiEXP::EvaluateFunction(double s) {
    return 1.0/(s-this->a);
}


#endif // METHODS_PSI_HPP_