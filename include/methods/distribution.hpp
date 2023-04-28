#ifndef METHODS_DISTRIBUTION_H_
#define METHODS_DISTRIBUTION_H_


#include <stdexcept>
#include <vector>
#include <map> 
#include <numeric>
#include <functional>
#include <random>
#include <fstream>
#include <algorithm>
#include <limits>
#include <iterator>

#include "exprtk.hpp"
#include "csv.h"

#include "io/echo.hpp"


class ExperimentalDistribution {
private:

    int                _seed;       ///< Seed
    std::random_device _device;     ///< Random device
    std::mt19937       _generator;  ///< Mersenne Twister PRNG

    std::uniform_real_distribution<double> _dist_uniform;  ///< Uniform distribution [0, 1]

    std::string                   _expression;  ///< Function expression: y = ax + b
    std::string                   _parameter;   ///< Function parameter:  x
    std::map<std::string, double> _constants;   ///< Function constants:  a, b

    int _num;  ///< Size of distribution

    std::vector<double> _x;    ///< x-axis values
    std::vector<double> _y;    ///< y-axis values
    std::vector<double> _pdf;  ///< Probability distribution function
    std::vector<double> _cdf;  ///< Cumulative distribution function

    double _x_min;    ///< Minimum _x value
    double _x_max;    ///< Maximum _x value
    double _y_min;    ///< Minimum _y value
    double _y_max;    ///< Maximum _y value
    double _pdf_min;  ///< Minimum _pdf_value
    double _pdf_max;  ///< Maximum _pdf_value
    double _mean;     ///< PDF expected value
    double _factor;   ///< Normalization factor

public:

    // Constructors

    ExperimentalDistribution(void);
    ExperimentalDistribution(std::vector<double> x, std::vector<double> y);
    ExperimentalDistribution(std::string expression, std::string parameter, std::map<std::string, double> constants, std::vector<double> x);

    // Getters

    const int GetSeed(void);

    const std::vector<double> GetX(void);
    const std::vector<double> GetY(void);
    const std::vector<double> GetPDF(void);
    const std::vector<double> GetCDF(void);

    const int GetSize(void);

    const double GetXMin(void);
    const double GetXMax(void);
    
    const double GetYMin(void);
    const double GetYMax(void);
    
    const double GetPDFMin(void);
    const double GetPDFMax(void);
    const double GetPDFMean(void);

    const double GetNormalizationFactor(void);

    const std::string GetExpression(void);
    const std::string GetParameter(void);
    const std::map<std::string, double> GetConstants(void);

    // Setters

    void SetSeed(int seed);

    void SetExpression(std::string expression);
    void SetConstants (std::map<std::string, double> constants);

    // Exceptions

    void VerifyShape(void);
    void VerifyMonotonicity(void);

    // Methods

    void CreateRandomGenerator(void);

    std::pair<double, double> CalculateMinMax(std::vector<double> values);
    double EstimateNumericalMax(std::vector<double> x_values, std::vector<double> y_values);

    void CalculateXMinMax(void);
    void CalculateYMinMax(void);
    void CalculatePDFMinMax(void);
    void CalculatePDFMean(void);
    void CalculateNormalizationFactor(void);

    void NormalizeFunction(void);
    void ParseExpression(void);

    // Samplers

    double SampleUniformDistribution(void);

    double SampleFunction(void);
    double SampleFunction(double z);

    // I/O

    void WriteFunction(const char *fname, const char delimiter=',');
    void ReadFunction(const char *fname);

};

class ExperimentalDistributionVector {
private:

    int _num;

public:

    ExperimentalDistributionVector(void);

};


// 
// Constructors
//

/// Constructs a new and empty ExperimentalDistribution object.
inline ExperimentalDistribution::ExperimentalDistribution(void) 
    : _seed(0), _num(0), _x(std::vector<double>(0)), _y(std::vector<double>(0)), _pdf(std::vector<double>(0)), _cdf(std::vector<double>(0)) {
}

/// Constructs a new ExperimentalDistribution object using an x and y vector. The x and y values will
/// be used to construct a pdf and cdf function.
inline ExperimentalDistribution::ExperimentalDistribution(std::vector<double> x, std::vector<double> y) 
    : _seed(0), _num(x.size()), _x(x), _y(y) {
    
    this->VerifyShape();
    this->VerifyMonotonicity();

    this->CalculateNormalizationFactor();
    this->NormalizeFunction();

    this->CalculateXMinMax();
    this->CalculateYMinMax();
    this->CalculatePDFMean();
}

/**
 * @brief Constructs a new ExperimentalDistribution object using an expression and its constants.
 * 
 * @param expression : the string containing the parsed function.
 * @param parameter  : the paramter symbol of the function.
 * @param constants  : the constants symbols and values.
 * @param x          : the vector with the paramters from which the expression will be evaluated.
 */
inline ExperimentalDistribution::ExperimentalDistribution(std::string expression, std::string parameter, std::map<std::string, double> constants, std::vector<double> x)
    : _seed(0), _expression(expression), _parameter(parameter), _constants(constants), _num(x.size()), _x(x) {

    this->VerifyMonotonicity();
    
    this->ParseExpression();
    this->CalculateNormalizationFactor();
    this->NormalizeFunction();

    this->CalculateXMinMax();
    this->CalculateYMinMax();
    this->CalculatePDFMean();
}

// 
// Getters
// 

inline const int ExperimentalDistribution::GetSeed(void) {
    return _seed;
}

inline const std::vector<double> ExperimentalDistribution::GetX(void) {
    return _x;
}

inline const std::vector<double> ExperimentalDistribution::GetY(void) {
    return _y;
}

inline const std::vector<double> ExperimentalDistribution::GetPDF(void) {
    return _pdf;
}

inline const std::vector<double> ExperimentalDistribution::GetCDF(void) {
    return _cdf;
}

inline const int ExperimentalDistribution::GetSize(void) {
    return _num;
}

inline const double ExperimentalDistribution::GetXMin(void) {
    return _x_min;
}

inline const double ExperimentalDistribution::GetXMax(void) {
    return _x_max;
}

inline const double ExperimentalDistribution::GetYMin(void) {
    return _y_min;
}

inline const double ExperimentalDistribution::GetYMax(void) {
    return _y_max;
}

inline const double ExperimentalDistribution::GetPDFMin(void) {
    return _pdf_min;
}

inline const double ExperimentalDistribution::GetPDFMax(void) {
    return _pdf_max;
}

inline const double ExperimentalDistribution::GetPDFMean(void) {
    return _mean;
}

inline const double ExperimentalDistribution::GetNormalizationFactor(void) {
    return _factor;
}

inline const std::string ExperimentalDistribution::GetExpression(void) {
    return _expression;
}

const std::string ExperimentalDistribution::GetParameter(void) {
    return _parameter;
}

const std::map<std::string, double> ExperimentalDistribution::GetConstants(void) {
    return _constants;
}


// 
// Setters
// 

inline void ExperimentalDistribution::SetSeed(int seed) {
    _seed = seed;
    this->CreateRandomGenerator();
}

inline void ExperimentalDistribution::SetExpression(std::string expression) {
    _expression = expression;
}

inline void ExperimentalDistribution::SetConstants(std::map<std::string, double> constants) {
    _constants = constants;
}

// 
// Methods
// 

inline void ExperimentalDistribution::VerifyShape(void) {
    if (_x.size() != _y.size())
        throw std::invalid_argument(
            "Mismatch between the shapes of x(" 
            + std::to_string(_x.size())  
            + ") and y(" 
            + std::to_string(_y.size()) 
            + ")"
        );

    if (_cdf.size() != 0)
        if ( (_cdf.size() !=_x.size()) || (_cdf.size() != _y.size()) )
            throw std::invalid_argument(
                "Mismatch between the shapes of x(" 
                + std::to_string(_x.size())  
                + "), y(" 
                + std::to_string(_y.size()) 
                + ") and y_cumsum("
                + std::to_string(_cdf.size()) 
                + ")"
            );
}

inline void ExperimentalDistribution::VerifyMonotonicity(void) {
    if (_x.size() != 0)
        if (!std::is_sorted(_x.begin(), _x.end()))
            throw std::invalid_argument("Invalid value, x has to be monotonic");
    
    if (_y.size() != 0)
        if (!std::is_sorted(_y.begin(), _y.end()))
            throw std::invalid_argument("Invalid value, y has to be monotonic");
    
    if (_cdf.size() != 0)
        if (!std::is_sorted(_cdf.begin(), _cdf.end()))
            throw std::invalid_argument("Invalid value, y_cumsum has to be monotonic");

}

inline void ExperimentalDistribution::CreateRandomGenerator(void) {
    _generator = std::mt19937(_device());
    if (_seed != 0) _generator.seed(_seed);
}

inline std::pair<double, double> ExperimentalDistribution::CalculateMinMax(std::vector<double> values) {
    std::pair<std::vector<double>::iterator, std::vector<double>::iterator> out 
        = std::minmax_element(values.begin(), values.end());

    return std::make_pair(*out.first, *out.second);
}

inline double ExperimentalDistribution::EstimateNumericalMax(std::vector<double> x_values, std::vector<double> y_values) {
    double eps  = std::numeric_limits<double>::epsilon();
    double step = 2.0;

    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double>   expression;
    exprtk::parser<double>       parser;
    
    double x;
    symbol_table.add_variable(_parameter, x);

    for (std::map<std::string, double>::iterator it = _constants.begin(); it != _constants.end(); ++it)
        symbol_table.add_variable(it->first, it->second);
        // std::cout << it->first << " => " << it->second << '\n';

    // for (std::pair<std::string, double> &var : _constants)
    //     symbol_table.add_variable(var.first, var.second);

    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);
    parser.compile(_expression, expression);

    double ta = 1.0;
    double tb = ta*step;

    x = ta;
    double pa = expression.value();
    x = tb;
    double pb = expression.value();

    double error = std::abs(pb-pa);

    while (error > eps) {
        tb *= step;
        x = tb;
        pb = expression.value();
        error = std::abs(pb-pa);
        ta = tb;
        pa = pb;
    }

    return ta;
}

inline void ExperimentalDistribution::CalculateXMinMax(void) {
    std::pair<double, double> val = this->CalculateMinMax(_x);
    _x_min = val.first;
    _x_max = val.second;
}

inline void ExperimentalDistribution::CalculateYMinMax(void) {
    std::pair<double, double> val = this->CalculateMinMax(_y);
    _y_min = val.first;
    _y_max = val.second;
}

inline void ExperimentalDistribution::CalculatePDFMinMax(void) {
    std::pair<double, double> val = this->CalculateMinMax(_pdf);
    _pdf_min = val.first;
    _pdf_max = val.second;
}

/**
 * @brief Calculates the PDF's expected value.
 * 
 * The expected value is calculated using the integral
 * $\operatorname{E} [X] = \int _{\mathbb {R}} x f(x) \, dx$
 * under our computational perspective.
 */
inline void ExperimentalDistribution::CalculatePDFMean(void) {
    _mean = 0.0;

    for (int i = 1; i < _num; ++i)
        _mean += 0.5*(_x[i]+_x[i-1])*0.5*(_pdf[i]+_pdf[i-1])*(_x[i]-_x[i-1]);
}

inline void ExperimentalDistribution::CalculateNormalizationFactor(void) {
    _factor = 0.0;

    for (int i = 1; i < _num; ++i)
        _factor += 0.5 * (_y[i]+_y[i-1]) * (_x[i]-_x[i-1]);

    _factor = 1.0 / _factor;
}

inline void ExperimentalDistribution::NormalizeFunction(void) {
    _cdf = std::vector<double>(_num, 0.0);
    _pdf = std::vector<double>(_num, 0.0);

    _pdf[0] = _factor * _y[0];
    _cdf[0] = _pdf[0] * (_x[1]-_x[0]);

    for (int i = 1; i < _num; ++i) {
        _pdf[i] = _y[i] * _factor;
        _cdf[i] = _cdf[i-1] + _pdf[i] * (_x[i]-_x[i-1]);
    }
}

inline void ExperimentalDistribution::ParseExpression(void) {
    // exp(-t/t2)*pow(1+t/t1,-1- beta)

    _y = std::vector<double>(_num);

    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double>   expression;
    exprtk::parser<double>       parser;
    
    double x;
    symbol_table.add_variable(_parameter, x);

    // for (std::pair<std::string, double> &var : _constants)
    //     symbol_table.add_variable(var.first, var.second);
    for (std::map<std::string, double>::iterator it = _constants.begin(); it != _constants.end(); ++it)
        symbol_table.add_variable(it->first, it->second);

    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);
    parser.compile(_expression, expression);

    for (int i=0; i<_num; ++i){
        x = _x[i];
        _y[i] = expression.value();
    }

}

// 
// Samplers
// 

/**
 * @brief Samples a random number from the interval [0, 1]. 
 * 
 * If the a seed was specified earlier, it's going to be used 
 * herein.
 * 
 * @return double : sample in the interval [0, 1]
 */
inline double ExperimentalDistribution::SampleUniformDistribution(void) {
    return _dist_uniform(_generator);
}

/**
 * @brief Samples a number using the CDF from the constructed 
 * experimental distribution.
 * 
 * @return double : the sampled value 
 */
inline double ExperimentalDistribution::SampleFunction(void) {
    return this->SampleFunction(this->SampleUniformDistribution());
}

/**
 * @brief Finds $x_i$, using the constructed CDF such that 
 * $\operatorname{CDF}(x_i) < z \leq \operatorname{CDF}(x_{i+1})$.
 * If the value lies in one of the $[x_i, x_{i+1}]$ extremes, 
 * its taken as it is, otherwise it's linearly interpolated from 
 * those values.
 * 
 * @param z       : random value from the interval [0, 1]
 * @return double : the sampled value
 */
inline double ExperimentalDistribution::SampleFunction(double z) {
    std::vector<double>::iterator lb = 
        std::lower_bound(_cdf.begin(), _cdf.end(), z);
    
    if (std::distance(_cdf.begin(), lb) == _num) 
        return _x[_num - 1];
    else
        return _x[lb - _cdf.begin()];

    // std::vector<double>::iterator ub = 
    //     std::upper_bound(_cdf.begin(), _cdf.end(), z);
    
    // double x0 = _x[ub - 1 - _cdf.begin()];
    // double x1 = _x[ub -     _cdf.begin()];

    // double y0 = _cdf[ub - 1 - _cdf.begin()];
    // double y1 = _cdf[ub -     _cdf.begin()];

    // return x0+(x1-x0)*(z-y0)/(y1-y0);
}

// 
// I/O
// 

inline void ExperimentalDistribution::WriteFunction(const char *fname, const char delimiter) {
    std::ofstream ofs(fname, std::ofstream::out);
    ofs << "x,y,pdf,cdf" << std::endl;
    
    for (int i = 0; i < _num; ++i) 
        ofs 
            << _x[i]   << delimiter
            << _y[i]   << delimiter
            << _pdf[i] << delimiter
            << _cdf[i] 
        << std::endl;

    ofs.close();
}

/**
 * @brief Reads distribution from file.
 * 
 * @param fname : input file name
 */
inline void ExperimentalDistribution::ReadFunction(const char *fname) {

    double x;
    double y;
    double cdf;

    io::CSVReader<3> in(fname);

    in.read_header(io::ignore_extra_column, "x", "y", "cdf");
    
    while(in.read_row(x, y, cdf)){
        _x.push_back(x);
        _y.push_back(y);
        _cdf.push_back(cdf);
    }

    this->VerifyMonotonicity();
}

// 
// Streams
// 

inline std::ostream &operator<<(std::ostream &os, const std::shared_ptr<ExperimentalDistribution> &f) {
    
    if (f->GetSeed() != 0)
        os << message("Seed",         f->GetSeed()               );
    
    os << message("Num. of samples",  f->GetSize()               );
    os << message("Min. x value",     f->GetXMin()               );
    os << message("Max. x value",     f->GetXMax()               );
    os << message("PDF mean value",   f->GetPDFMean()            );
    os << message("PDF norm. factor", f->GetNormalizationFactor());

    if (f->GetExpression().size() > 0) {
        os << message("Expression",  f->GetExpression()         );
        os << message("Parameter",   f->GetParameter()          );
    }

    for (std::pair<std::string, double> var : f->GetConstants())
        os << message("Const. " + var.first,     var.second     );

    return os;
}


// 
// Methods
// 

inline std::vector<std::shared_ptr<ExperimentalDistribution>> CreateLocalDistribution(std::string expression, std::string parameter, std::vector<std::map<std::string, double>> constants, std::vector<double> x) {

    std::vector<std::shared_ptr<ExperimentalDistribution>> out(constants.size());

    for (int i=0; i<constants.size(); ++i) {
        out[i] = std::make_shared<ExperimentalDistribution>(expression, parameter, constants[i], x);
    }

    return std::move(out);
}

#endif  // METHODS_DISTRIBUTION_H_
