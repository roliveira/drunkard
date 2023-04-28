#ifndef PARTICLE_PARTICLE_H_
#define PARTICLE_PARTICLE_H_


#include <vector>
#include <random>
#include <memory>
#include <numeric>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "geometry/geometry.hpp"
#include "utils/index.hpp"
#include "utils/rng.hpp"

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;


class Particle
{

public:

    bool _move_list = false;

    int _id;

    int _i, _j, _k;
    int _i0, _j0, _k0;

    int _dir;

    float _x, _y, _z;
    float _x0, _y0, _z0;

    xt::xarray<double> xg;
    xt::xarray<float> dx = {0.0, 0.0, 0.0};

    int _step;
    float _t0, _t, _tt;
    float _taut, _taur;
    float _dtt, _dtr;

    double v = 0.0;

    double mols = 0.0;

    Reaktoro::VectorXd amounts;

    bool _valid;
    int _type;

    ~Particle(void);

    Particle(void);
    Particle(int id);
    Particle(int id, float x, float y, float z);

    void SetX0(float x0, float y0, float z0);
    void SetI0(int i0, int j0, int k0);
    void SetT0(float t0, float tt0);

    void SetX(float x, float y, float z);
    void SetI(int i, int j, int k);

    xt::xarray<int> GetSourceI(void);
    xt::xarray<int> GetTargetI(void);

};

Particle::~Particle(void)
{}

Particle::Particle(void)
{}

Particle::Particle(int id): 
_id(id)
{}

Particle::Particle(int id, float x, float y, float z): 
_id(id)
{
    this->_x = x;
    this->_y = y;
    this->_z = z;
    
    this->_x0 = x;
    this->_y0 = y;
    this->_z0 = z;
}

void Particle::SetX0(float x0, float y0, float z0)
{
    this->_x0 = x0;
    this->_y0 = y0;
    this->_z0 = z0;

    this->xg = xt::xarray<float>{x0, y0, z0};
}

void Particle::SetI0(int i0, int j0, int k0)
{   
    this->_i0 = i0;
    this->_j0 = j0;
    this->_k0 = k0;
}

void Particle::SetT0(float t0, float tt0)
{
    this->_t0  = t0;
    this->_t   = t0;
    this->_tt  = tt0;
    this->_dtt = tt0;
}

void Particle::SetX(float x, float y, float z)
{
    this->_x = x;
    this->_y = y;
    this->_z = z;
}

void Particle::SetI(int i, int j, int k)
{
    this->_i = i;
    this->_j = j;
    this->_k = k;
}

xt::xarray<int> Particle::GetSourceI(void)
{
    return xt::xarray<int>({ this->_i, this->_j, this->_k });
}

xt::xarray<int> Particle::GetTargetI(void)
{
    xt::xarray<int> ind { this->_i, this->_j, this->_k };

    xt::xarray<int> vector_dir = 
        xt::view(direction_vector(), this->_dir, xt::all());

    return ind + vector_dir;
}


#endif  // PARTICLE_PARTICLE_H_
