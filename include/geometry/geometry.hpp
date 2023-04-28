#ifndef GEOMETRY_GEOMETRY_H_
#define GEOMETRY_GEOMETRY_H_


#include <iostream>
#include <memory>
#include <vector>
#include <array>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <stdexcept>
#include <random>

#include "config.hpp"
#include "io/echo.hpp"
#include "io/reader.hpp"
#include "io/writer.hpp"
#include "utils/index.hpp"
#include "utils/operations.hpp"

#include "xtensor/xarray.hpp"
#include "xtensor/xbuilder.hpp"


class Geometry {

public:

    int _length;
    std::vector<int> _shape;
    std::vector<double> _size;
    std::vector<double> _origin;
    std::vector<double> _spacing;
    
    xt::xarray<double> spacing;

    // Constructors

    Geometry(void);
    Geometry(std::vector<int> shape);
    Geometry(std::vector<double> size );
    Geometry(std::vector<int> shape, std::vector<double> size);
    Geometry(std::vector<int> shape, std::vector<double> size, std::vector<double> origin);

    // Methods
    
    void CalculateLength(void);
    void CalculateSpacing(void);
    void AddGhostLayer(void);
    void CalculateLocalPeclet(void);

    std::vector<double> CreateLinearField(double _x0_val, double _x1_val);
    std::vector<double> DisturbField(std::vector<double> _field, double _mag);

    std::vector<int> Ind2Sub(int ind);
    int Sub2Ind(std::vector<int> sub);

    double CalculateDistance(xt::xarray<int> i, xt::xarray<int> j);

    // Getters

    double GetDistance(int i, int j);
    double GetDistance(std::vector<int> i, std::vector<int> j);

    std::vector<double> GetLocation(int linear_index );
    std::vector<double> GetLocation(int i, int j, int k );
    std::vector<double> GetLocation(std::vector<int> grid_subindex );
    xt::xarray <double> GetLocation(xt::xarray<int> subind);

    int GetReflectedGridIndex(            int  _linear_index , int _lag=1);
    int GetReflectedGridIndex(std::vector<int> _grid_subindex, int _lag=1);

    std::vector<int> GetReflectedGridSubIndex(            int  _linear_index , int _lag=1);
    std::vector<int> GetReflectedGridSubIndex(std::vector<int> _grid_subindex, int _lag=1);

    xt::xarray<int> GetBouncedGridSubindex(xt::xarray<int> subind);
    xt::xarray<int> GetReflectedGridSubIndex(xt::xarray<int> subind);

    std::vector<int> GetFaceNeighboursIndices(            int  _linear_index );
    std::vector<int> GetFaceNeighboursIndices(std::vector<int> _grid_subindex);

    // Comparisons

    bool IsBoundary(int ind);
    bool IsBoundary(int i, int j, int k);
    bool IsBoundary(std::vector<int> subind);
    bool IsBoundary(xt::xarray<int> subind);

    bool IsXm(xt::xarray<int> subind);
    bool IsXp(xt::xarray<int> subind);
    bool IsYm(xt::xarray<int> subind);
    bool IsYp(xt::xarray<int> subind);
    bool IsZm(xt::xarray<int> subind);
    bool IsZp(xt::xarray<int> subind);

};

//
// Constructors
//

Geometry::Geometry(void)
: _length(0), _shape({}), _size({}), _origin({}), _spacing({}) 
{}

Geometry::Geometry(std::vector<int> shape)
: Geometry(shape, std::vector<double>(shape.size(), 1.0), std::vector<double>(shape.size(), 0.0f))
{}

Geometry::Geometry(std::vector<double> size)
: Geometry(std::vector<int>(size.size(), 1), size, std::vector<double>(size.size(), 0.0f))
{}

Geometry::Geometry(std::vector<int> shape, std::vector<double> size)
: Geometry(shape, size, std::vector<double>(shape.size(), 0.0f))
{}

Geometry::Geometry(std::vector<int> shape, std::vector<double> size, std::vector<double> origin) 
: _shape(shape), _size(size), _origin(origin) 
{
    this->CalculateSpacing();
    this->AddGhostLayer();
    this->CalculateLength();
}

//
// Methods
//

void Geometry::AddGhostLayer(void) 
{
    std::transform(
        this->_shape.begin(), 
        this->_shape.end(), 
        this->_shape.begin(), 
        bind2nd(std::plus<int>(), 2)
    );   
}

void Geometry::CalculateLength(void) 
{
    this->_length = std::accumulate(
        this->_shape.begin(), 
        this->_shape.end(), 
        1, 
        std::multiplies<int>()
    );
}

void Geometry::CalculateSpacing(void) 
{
    this->_spacing = std::vector<double>(this->_size);
    
    for (int i=0; i<this->_spacing.size(); ++i)
    {
        this->_spacing[i] /= (this->_shape[i]-1);
    }

    this->spacing = xt::adapt(this->_spacing, {3});
}

std::vector<double> Geometry::CreateLinearField(double x0_val, double x1_val) 
{
    std::vector<double> out(this->_length, x0_val);
    double a = (x0_val-x1_val)/(0.0-this->_shape[0]);

    for (int i = 0; i < out.size(); ++i)
        out[i] += a*(ind2sub(i, this->_shape)[0]-x1_val);

    return out;
}

std::vector<double> Geometry::DisturbField(std::vector<double> _field, double _mag) 
{
    std::random_device device;
    std::mt19937 generator(device());

    std::uniform_real_distribution<> dist = 
        std::uniform_real_distribution<>(-1, 1);

    for(double &i : _field)
        i += dist(generator) * _mag;

    return _field;
}

int Geometry::Sub2Ind(std::vector<int> sub) 
{
    return sub2ind(sub, this->_shape);
}

std::vector<int> Geometry::Ind2Sub(int ind) 
{
    return ind2sub(ind, this->_shape);
}

double Geometry::CalculateDistance(xt::xarray<int> i, xt::xarray<int> j)
{
    return VectorMagnitude(this->GetLocation(i), this->GetLocation(j));
}


//
// Setters
//

//
// Getters
//


double Geometry::GetDistance(int i, int j) 
{
    return VectorMagnitude(this->GetLocation(i), this->GetLocation(j));
}

double Geometry::GetDistance(std::vector<int> i, std::vector<int> j)
{
    return VectorMagnitude(this->GetLocation(i), this->GetLocation(j));
}

// Get Location

std::vector<double> Geometry::GetLocation(int ind) 
{
    return this->GetLocation(ind2sub(ind, this->_shape));
}

std::vector<double> Geometry::GetLocation(int i, int j, int k) 
{
    return this->GetLocation(std::vector<int>({i, j, k}));
}

std::vector<double> Geometry::GetLocation(std::vector<int> subind) 
{
    std::vector<double> location(3);
    
    for (int i = 0; i < 3; ++i) 
        location[i] = this->_spacing.at(i)*static_cast<double>(subind.at(i));

    return location;
}

xt::xarray<double> Geometry::GetLocation(xt::xarray<int> subind)
{
    xt::xarray<double> location = xt::zeros<double>({3});
    
    for (int i = 0; i < 3; ++i) 
        location[i] = this->_spacing.at(i)*static_cast<double>(subind.at(i));

    return location;
}

// Get Reflected Index

int Geometry::GetReflectedGridIndex(int ind, int lag) 
{
    return mirror_index(ind, this->_shape, lag);
}

int Geometry::GetReflectedGridIndex(std::vector<int> subind, int lag) 
{
    return this->GetReflectedGridIndex(sub2ind(subind, this->_shape), lag);
}

// Get Reflected Sub-index

std::vector<int> Geometry::GetReflectedGridSubIndex(int ind, int lag) 
{
    return this->GetReflectedGridSubIndex(ind2sub(ind, this->_shape), lag);
}

std::vector<int> Geometry::GetReflectedGridSubIndex(std::vector<int> subind, int lag) 
{
    return mirror_subindex(subind, this->_shape, lag);
}

xt::xarray<int> Geometry::GetReflectedGridSubIndex(xt::xarray<int> subind)
{
    xt::xarray<int> out = xt::zeros<int>({3});

    for (int i=0; i<3; ++i) {
        if (subind(i) <= 1                ) out(i) = this->_shape[i]-3;
        if (subind(i) >= this->_shape[i]-2) out(i) = 2;
    }

    return out;
}

xt::xarray<int> Geometry::GetBouncedGridSubindex(xt::xarray<int> subind)
{
    xt::xarray<int> out = subind;

    for (int i=0; i<3; ++i) {
        if (subind(i) <= 1                ) out(i) = 2;
        if (subind(i) >= this->_shape[i]-2) out(i) = this->_shape[i]-3;
    }

    return out;
}

// Neighbours

std::vector<int> Geometry::GetFaceNeighboursIndices(int ind) 
{
    return this->GetFaceNeighboursIndices(ind2sub(ind, this->_shape));
}

std::vector<int> Geometry::GetFaceNeighboursIndices(std::vector<int> subind) 
{
    return face_neighbours_index(3, sub2ind(subind, this->_shape), this->_shape);
}


// Comparisons

bool Geometry::IsBoundary(int ind) 
{
    return this->IsBoundary(ind2sub(ind, this->_shape));
}

bool Geometry::IsBoundary(int i, int j, int k) 
{
    return is_boundary(std::vector<int>({i, j, k}), this->_shape);
}

bool Geometry::IsBoundary(std::vector<int> subind) 
{
    return is_boundary(subind, this->_shape);
}

bool Geometry::IsBoundary(xt::xarray<int> subind) 
{
    return this->IsBoundary(subind(0), subind(1), subind(2));
}


bool Geometry::IsXm(xt::xarray<int> subind)
{
    return subind(0) <= 1;
}

bool Geometry::IsXp(xt::xarray<int> subind)
{
    return subind(0) >= this->_shape[0]-2;
}

bool Geometry::IsYm(xt::xarray<int> subind)
{
    return subind(1) <= 1;
}

bool Geometry::IsYp(xt::xarray<int> subind)
{
    return subind(1) >= this->_shape[1]-2;
}

bool Geometry::IsZm(xt::xarray<int> subind)
{
    return subind(1) <= 1;
}

bool Geometry::IsZp(xt::xarray<int> subind)
{
    return subind(2) >= this->_shape[2]-2;
}


//
// I/O
//

// Streams 

std::ostream &operator<<(std::ostream &os, std::shared_ptr<Geometry> const &gs)
{
    os << message("Shape",       gs->_shape);
    os << message("Origin [m]",  gs->_origin);
    os << message("Size [m]",    gs->_size);
    os << message("Spacing [m]", gs->_spacing);
    
    return os;
}


#endif  // GEOMETRY_GEOMETRY_H_
