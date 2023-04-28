#ifndef UTILS_RNG_HPP_
#define UTILS_RNG_HPP_


#include <random>


class RandomNumberGenerator
{

public:

    int                _seed;
    std::random_device _device;
    std::mt19937       _generator;

    std::uniform_real_distribution<double> _uniform;

    RandomNumberGenerator(void);
    RandomNumberGenerator(int seed);

    float Sample(void);

};


inline RandomNumberGenerator::RandomNumberGenerator(void)
{ 
    this->_generator = std::mt19937(this->_device());
    this->_uniform   = std::uniform_real_distribution<double>(0.0, 1.0);
}

inline RandomNumberGenerator::RandomNumberGenerator(int seed)
: _seed(seed)
{
    RandomNumberGenerator();
    this->_generator.seed(this->_seed);
}

inline float RandomNumberGenerator::Sample(void)
{
    return this->_uniform(this->_generator);
}


typedef RandomNumberGenerator RNG;


#endif // UTILS_RNG_HPP_
