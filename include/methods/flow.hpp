
#ifndef METHODS_FLOW_H_
#define METHODS_FLOW_H_


#include <iostream>
#include <vector>
#include <iterator>

#include "io/writer.hpp"
#include "geometry/geometry.hpp"


enum PRESSURE_BC {
    NONE = 0,
    PIN  = 1,
    POUT = 2
};


class FlowMethod : public Geometry {
public:

    std::vector<PRESSURE_BC> _flag_bc;

    // Constructor

    FlowMethod(void);
    FlowMethod(std::vector<int> shape);
    FlowMethod(std::vector<double> size );
    FlowMethod(std::vector<int> shape, std::vector<double> size);
    FlowMethod(std::vector<int> shape, std::vector<double> size, std::vector<double> origin);

    // Setters

    void SetBoundaryCondition(int ind, PRESSURE_BC flag_bc);

    // Methods

    std::vector<double> CreateCorrelatedField(double mean, double stddev, int num, int seed);
    std::vector<double> CreateNormalField    (double mean, double stddev, int num, int seed);
    std::vector<double> CreateLogNormalField (double mean, double stddev, int num, int seed);

    // I/O

    void WriteScalarField(const char *_fname, std::vector<double> _data);
    void WriteVectorField(const char *_fname, std::vector<double> _data);
    void WriteTensorField(const char *_fname, std::vector<double> _data);

    std::vector<double> ReadScalarField(const char *_fname);
    std::vector<double> ReadVectorField(const char *_fname);
    std::vector<double> ReadTensorField(const char *_fname);
};


// Constructors

inline FlowMethod::FlowMethod(void)
: Geometry() 
{}

inline FlowMethod::FlowMethod(std::vector<int> shape)
: Geometry(shape) 
{}

inline FlowMethod::FlowMethod(std::vector<double> size) 
: Geometry(size) 
{}

inline FlowMethod::FlowMethod(std::vector<int> shape, std::vector<double> size)
: Geometry(shape, size) 
{}

inline FlowMethod::FlowMethod(std::vector<int> shape, std::vector<double> size, std::vector<double> origin) 
: Geometry(shape, size, origin) 
{
    _flag_bc = std::vector<PRESSURE_BC>(this->_length, NONE);
}

// Setters

inline void FlowMethod::SetBoundaryCondition(int ind, PRESSURE_BC flag_bc) 
{
    _flag_bc[ind] = flag_bc;
}

// Methods

inline std::vector<double> FlowMethod::CreateCorrelatedField(double val_mean, double val_stddev, int num_val, int seed) 
{
    int lag = 3;

    std::random_device device;
    std::mt19937 generator = std::mt19937(device());
    if (seed != 0) generator.seed(seed);

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::vector<double> out = 
        this->CreateNormalField(val_mean, val_stddev, num_val, seed);

    for (int iz=0; iz<this->_shape.at(2); ++iz) 
    {
        int iiz = iz;

        for (int iy=0; iy<this->_shape.at(1); ++iy) 
        {
            int iiy = iy;

            for (int ix=0; ix<this->_shape.at(0); ++ix) 
            {
                int iix = ix;

                int ind = this->Sub2Ind({ix, iy, iz});

                if (ix<lag || ix>this->_shape.at(0)-1-lag ||
                    iy<lag || iy>this->_shape.at(1)-1-lag || 
                    iz<lag || iz>this->_shape.at(2)-1-lag) 
                    {
                        std::vector<int> rind = this->GetReflectedGridSubIndex({iix, iiy, iiz}, lag);
                        iix = rind[0];
                        iiy = rind[1];
                        iiz = rind[2];
                    }

                for (int i=-lag; i<=lag; ++i)
                {
                    int iind = this->Sub2Ind({iix+i, iiy, iiz});
                    out[iind] += val_stddev * dist(generator) / (6*(lag+1));

                    iind = this->Sub2Ind({iix, iiy+i, iiz});
                    out[iind] += val_stddev * dist(generator) / (6*(lag+1));

                    iind = this->Sub2Ind({iix, iiy, iiz+i});
                    out[iind] += val_stddev * dist(generator) / (6*(lag+1));
                }
            }
        }
    }

    return out;
}


inline std::vector<double> FlowMethod::CreateNormalField(double val_mean, double val_stddev, int num_val, int seed) 
{
    std::random_device device;
    std::mt19937 generator = std::mt19937(device());

    if (seed != 0) generator.seed(seed);

    std::normal_distribution<double> dist{val_mean, val_stddev};

    std::vector<double> out(num_val);
    std::generate(
        out.begin(), 
        out.end(), 
        [&]() { 
            double z;
            do {z = dist(generator);} while (z<0);
            return z;
        }
    );

    return out;
}

inline std::vector<double> FlowMethod::CreateLogNormalField(double val_mean, double val_stddev, int num_val, int seed) 
{
    std::random_device device;
    std::mt19937 generator = std::mt19937(device());

    if (seed != 0) 
        generator.seed(seed);

    std::lognormal_distribution<double> dist{val_mean, val_stddev};

    std::vector<double> out(num_val);
    std::generate(out.begin(), out.end(), [&]() { return dist(generator); });

    return out;
}

// I/O

inline void FlowMethod::WriteScalarField(const char *_fname, std::vector<double> _data) {
    ToVTI(_fname, _data, this->_shape, this->_spacing, this->_origin);
}

inline void FlowMethod::WriteVectorField(const char *_fname, std::vector<double> _data) {
    Vector2VTI(_fname, _data, this->_shape, this->_spacing, this->_origin);
}

inline void FlowMethod::WriteTensorField(const char *_fname, std::vector<double> _data) {
    Tensor2VTI(_fname, _data, this->_shape, this->_spacing, this->_origin);
}

#endif  // METHODS_FLOW_H_
