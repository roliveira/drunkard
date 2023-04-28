
#ifndef METHODS_CONTINUOUS_HPP_
#define METHODS_CONTINUOUS_HPP_


#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <limits>
#include <limits>

#include "histogram.hpp"

#include "methods/transport.hpp"
#include "methods/distribution.hpp"
#include "methods/psi.hpp"
#include "utils/root.hpp"
#include "utils/interpolation.hpp"
#include "utils/quadrature.hpp"
#include "utils/hist.hpp"
#include "io/writer.hpp"


class ContinuousTimeRandomWalk : public TransportMethod {
private:

    // Distributions

    std::vector<std::shared_ptr<ExperimentalDistribution>> dist_sampling;
    std::vector<Psi> psi_dist;

    // Direction

    int                 ndirections;
    std::vector<int   > direction_vector;
    std::vector<double> direction_prob;
    std::vector<double> direction_prob_cumm;

public:

    // Constructors

    ContinuousTimeRandomWalk(void);
    ContinuousTimeRandomWalk(int _seed);
    ContinuousTimeRandomWalk(double _time_initial, double _time_delta, double _time_max, std::string _fprefix);
    ContinuousTimeRandomWalk(double _time_initial, double _time_delta, double _time_max, std::string _fprefix, int _seed);

    // Setters

    void SetSamplingDistribution(std::shared_ptr<ExperimentalDistribution> f);
    void SetSamplingDistribution(std::vector<std::shared_ptr<ExperimentalDistribution>> f);
    void SetSamplingDistribution(std::vector<Psi> f);

    void SetInitialCondition(void);
    void SetNextLocation(int particle_index, std::vector<int> next_location);

    // Methods

    // std::vector<double> CreateVelocityFieldPerDomain(std::vector<double> v, std::vector<int> dir);
    
    void CreateLocationFunction(double _t1, double _t2);
    void CreateLocationFunction(std::vector<double> _velocity, double _t2);
    // void CreateLocationFunctionPerDomain(std::vector<double> _t1, std::vector<double> _t2);

    void IterateParticles(void);
    void IterateParticlesSingleStep(void);

    // Samplers

    double SamplePsiFunction(int _particle_index);
    std::vector<int> SampleNextLocation(int particle_index);

    // I/O

    void WriteLocationTensor(void);

    void WriteSampledPsiFunction     (int n);
    void WriteSampledLocationFunction(int n);
    void WriteDisplacementProbability(const char *fname);

};

//
// Constructor
// 

inline ContinuousTimeRandomWalk::ContinuousTimeRandomWalk(void) 
: TransportMethod()
{}

inline ContinuousTimeRandomWalk::ContinuousTimeRandomWalk(int _seed) 
: TransportMethod(_seed)
{}

inline ContinuousTimeRandomWalk::ContinuousTimeRandomWalk(double _time_initial, double _time_delta, double _time_max, std::string _fprefix)
: TransportMethod(_time_initial, _time_delta, _time_max, _fprefix)
{}

inline ContinuousTimeRandomWalk::ContinuousTimeRandomWalk(double _time_initial, double _time_delta, double _time_max, std::string _fprefix, int _seed)
: TransportMethod(_seed, _time_initial, _time_delta, _time_max, _fprefix)
{}

//
// Setters
//

inline void ContinuousTimeRandomWalk::SetSamplingDistribution(std::shared_ptr<ExperimentalDistribution> f) {
    dist_sampling = { std::move(f) };
}

inline void ContinuousTimeRandomWalk::SetSamplingDistribution(std::vector<std::shared_ptr<ExperimentalDistribution>> f) {
    for (std::shared_ptr<ExperimentalDistribution> &v : f)
        dist_sampling.push_back(std::move(v));
}

inline void ContinuousTimeRandomWalk::SetSamplingDistribution(std::vector<Psi> f)
{
    psi_dist = f;
}


inline void ContinuousTimeRandomWalk::SetInitialCondition(void) 
{
    this->GetParticleSystem()->SetDeltaTime(this->GetDeltaTime());
    this->GetParticleSystem()->SetInitialStep(0);
    this->GetParticleSystem()->SetInitialTime(this->GetInitialTime());

    for (int ip = 0; ip < this->GetParticleSystem()->GetNumberOfParticles(); ++ip) 
    {
        this->GetParticleSystem()->SetTransitTime(ip, this->SamplePsiFunction(ip));
    }

    this->GetParticleSystem()->SetInitialVariance();
    this->GetParticleSystem()->CalculateParticleSummary();

    old_conc = std::vector<double>(this->GetGeometrySystem()->GetLinearShape(), 0.0);
}

inline void ContinuousTimeRandomWalk::SetNextLocation(int particle_index, std::vector<int> next_direction) {
    std::vector<int>    subind          = this->GetParticleSystem()->GetGridSubIndex(particle_index);
    std::vector<double> location        = this->GetGeometrySystem()->GetLocation    (subind        );
    std::vector<double> location_global = this->GetParticleSystem()->GetPosition    (particle_index);

    // Move particle to the next direction

    for (int i = 0; i < next_direction.size(); ++i)
        location_global[i] = 
            location_global[i] + next_direction[i] * this->GetGeometrySystem()->GetSpacing(i);

    // Increment the sub-indeces 

    for (int i = 0; i < next_direction.size(); ++i)
        subind[i] = subind[i] + next_direction[i];

    // Correct the locations and indices for the particles in the ghost layer

    if (this->GetGeometrySystem()->IsBoundary(subind)) {
        subind   = this->GetGeometrySystem()->GetReflectedGridSubIndex(subind);
        location = this->GetGeometrySystem()->GetLocation(subind);
        
        /*
        if (this->GetBoundaryCondition() == OUT) {
            this->GetParticleSystem()->SetParticleState(particle_index, STOP);
        }

        if (this->GetBoundaryCondition() == REV) {
            subind   = this->GetGeometrySystem()->GetReflectedGridSubIndex(subind);
            location = this->GetGeometrySystem()->GetLocation(subind);
        }

        if (this->GetBoundaryCondition() == INJ) {
            subind   = this->GetGeometrySystem()->GetReflectedGridSubIndex(subind);
            location = this->GetGeometrySystem()->GetLocation(subind);
            location_global = location;
        }
        */
    }

    this->GetParticleSystem()->SetGridSubIndex   ( particle_index, subind          ); 
    this->GetParticleSystem()->SetPositionLocal  ( particle_index, location        );
    this->GetParticleSystem()->SetPositionGlobal ( particle_index, location_global );
}

//
// Getters
//

// 
// Methods
//

inline void ContinuousTimeRandomWalk::CreateLocationFunction(double _t1, double _t2) {
    
    int dimension = this->GetDimension();
    double pe = 2 * _t2 / _t1;

    direction_vector  = face_vectors(dimension);
    ndirections       = 2 * dimension;
    direction_prob    = std::vector<double>(ndirections, 1.0);
    direction_prob[0] = pe / (std::exp(+pe) - 1);
    direction_prob[1] = pe / (1 - std::exp(-pe));

    direction_prob_cumm    = std::vector<double>(ndirections, 0.0);
    direction_prob_cumm[0] = direction_prob[0];

    for (int i = 1; i < direction_prob_cumm.size(); ++i) 
        direction_prob_cumm[i] = direction_prob[i] + direction_prob_cumm[i-1];

    for (double &i : direction_prob_cumm) 
        i /= direction_prob_cumm[direction_prob_cumm.size()-1];

}

inline void ContinuousTimeRandomWalk::CreateLocationFunction(std::vector<double> _velocity, double _t2) {
    
    float eps     = std::numeric_limits<float>::epsilon();
    int dimension = this->GetDimension();
    ndirections   = 2 * dimension;
    int nsize     = _velocity.size() / ndirections;
    
    direction_vector    = face_vectors(dimension);
    direction_prob      = std::vector<double>(_velocity.size(), 0.0);
    direction_prob_cumm = std::vector<double>(direction_prob.size(), 0.0);

    for (int i = 0; i < nsize; ++i) {
        if (this->GetGeometrySystem()->IsBoundary(i)) continue;

        for (int j=0; j<ndirections; ++j) {
            double v  =  _velocity[i*ndirections + j];
            double l  = this->GetGeometrySystem()->GetSize(j / 2);
            double t1 = std::abs(l / v);
            double pe = 2 * _t2 / t1;

            // pe ~~ 0 : pe / (1 - e^-pe) -> 1,
            //           pe / (e^+pe - 1) -> 1.
            if (std::abs(pe) < eps * 1.0){
                direction_prob[i*ndirections + j] = 1.0;
            }
            // pe >> 1 : pe / (1 - e^-pe) -> pe,
            //           pe / (e^+pe - 1) -> 0.
            else if (pe > 50) {
                if (v > 0) direction_prob[i*ndirections + j] = pe;
                else       direction_prob[i*ndirections + j] = 0.0;
            }
            else {
                if (v > 0) direction_prob[i*ndirections + j] = pe / (1 - std::exp(-pe));
                else       direction_prob[i*ndirections + j] = pe / (std::exp(+pe) - 1);
            }
        }
    }

    for (int i = 0; i < nsize; ++i) {
        if (this->GetGeometrySystem()->IsBoundary(i)) continue;

        std::partial_sum(
            direction_prob.begin()       +  i      * ndirections, 
            direction_prob.begin()       + (i + 1) * ndirections, 
            direction_prob_cumm.begin()  +  i      * ndirections, 
            std::plus<double>()
        );

        for (int j = 0; j < ndirections; ++j) {
            direction_prob[i*ndirections + j] /=
                direction_prob_cumm[(i + 1)*ndirections - 1];
                
            direction_prob_cumm[i*ndirections + j] /= 
                direction_prob_cumm[(i + 1)*ndirections - 1];
        }
    }

    for (int i = 0; i < nsize; ++i) {
        if (!this->GetGeometrySystem()->IsBoundary(i)) continue;

        for (int j = 0; j < ndirections; ++j) {
            int ir = this->GetGeometrySystem()->GetReflectedGridIndex(i);

            direction_prob[i*ndirections + j] = direction_prob[ir*ndirections + j];
            direction_prob_cumm[i*ndirections + j] = direction_prob_cumm[ir*ndirections + j];
        }
    }

}

inline void ContinuousTimeRandomWalk::IterateParticles(void) {
    for (int ip = 0; ip < this->GetParticleSystem()->GetNumberOfParticles(); ++ip) {
        if (this->GetParticleSystem()->GetState(ip) == 0) continue;
        
        while (this->GetParticleSystem()->GetPhysicalTime(ip) < this->GetCurrentTime()) {
            if (this->GetParticleSystem()->GetPhysicalTime(ip) + this->GetParticleSystem()->GetTransitTime(ip) > this->GetCurrentTime() + this->GetDeltaTime()) {
                this->GetParticleSystem()->SetPhysicalTimeIncrement(ip, this->GetDeltaTime(), -this->GetDeltaTime());
            }
            else {
                this->GetParticleSystem()->SetPhysicalTimeIncrement(ip, this->GetParticleSystem()->GetTransitTime(ip), 0.0);
                this->GetParticleSystem()->SetTransitTime(ip, this->SamplePsiFunction(ip));
                this->SetNextLocation(ip, this->SampleNextLocation(ip));
            }
        }
    }

    this->GetParticleSystem()->IncrementStep();
    this->GetParticleSystem()->IncrementTime(this->GetDeltaTime());
    this->GetParticleSystem()->CalculateParticleSummary();

    this->SetCurrentTime(this->GetCurrentTime() + this->GetDeltaTime());
}

inline void ContinuousTimeRandomWalk::IterateParticlesSingleStep(void) 
{
    double dt;

    for (int ip = 0; ip < this->GetParticleSystem()->GetNumberOfParticles(); ++ip) 
    { 
        dt = this->GetParticleSystem()->GetTransitTime(ip);

        this->GetParticleSystem()->SetPhysicalTimeIncrement(ip, dt, 0.0);
        this->GetParticleSystem()->SetTransitTime(ip, this->SamplePsiFunction(ip));
        this->SetNextLocation(ip, this->SampleNextLocation(ip));
    }

    this->GetParticleSystem()->IncrementStep();
    this->GetParticleSystem()->IncrementTime(dt);
    this->GetParticleSystem()->CalculateParticleSummary();

    this->SetCurrentTime(this->GetCurrentTime() + dt);
}

//
// Samplers
//

inline double ContinuousTimeRandomWalk::SamplePsiFunction(int _particle_index) 
{
    int ind = this->GetGeometrySystem()->Sub2Ind(
        this->GetParticleSystem()->GetGridSubIndex(_particle_index)
    );

    // return dist_sampling[ind]->SampleFunction();
    return psi_dist[ind].SampleTime(this->Sample());
}

inline std::vector<int> ContinuousTimeRandomWalk::SampleNextLocation(int _particle_index)
{
    int ind = this->GetGeometrySystem()->Sub2Ind(
        this->GetParticleSystem()->GetGridSubIndex(_particle_index)
    );

    std::vector<double> temp_direction_prob_cumm = 
        std::vector<double>(
            direction_prob_cumm.begin()+ ind   *ndirections,
            direction_prob_cumm.begin()+(ind+1)*ndirections
        );

    std::vector<double>::iterator lb = 
        std::lower_bound(
            temp_direction_prob_cumm.begin(),
            temp_direction_prob_cumm.end(),
            this->Sample()
        );

    int dimension = this->GetDimension();
    int idx = std::distance(temp_direction_prob_cumm.begin(), lb);

    if (idx == ndirections)
    {
        return std::vector<int>(
            &direction_vector[dimension* ndirections-1   ],
            &direction_vector[dimension*(ndirections-1+1)]
        );
    }
    else 
    {
        return std::vector<int>(
            &direction_vector[dimension* idx   ],
            &direction_vector[dimension*(idx+1)]
        );
    }
}

//
// I/O
//

inline void ContinuousTimeRandomWalk::WriteLocationTensor(void) {
    Tensor2VTI(
        (this->GetFilePrefix() + "_loc_prob.vti").c_str(),
        direction_prob,
        this->GetGeometrySystem()->GetShape(),
        this->GetGeometrySystem()->GetSpacing(),
        this->GetGeometrySystem()->GetOrigin()
    );

    Tensor2VTI(
        (this->GetFilePrefix() + "_loc_prob_cumm.vti").c_str(),
        direction_prob_cumm,
        this->GetGeometrySystem()->GetShape(),
        this->GetGeometrySystem()->GetSpacing(),
        this->GetGeometrySystem()->GetOrigin()
    );
}

// Streams

inline std::ostream &operator<<(std::ostream &os, const std::shared_ptr<ContinuousTimeRandomWalk> &tm) 
{
    os << message("Transport method",       "'Continuous Time Random Walk'");
    os << static_cast<const std::shared_ptr<TransportMethod>>(tm);
    return os;
}


#endif  // METHODS_CONTINUOUS_HPP_
