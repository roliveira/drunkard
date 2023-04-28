
#ifndef METHODS_TRANSPORT_H_
#define METHODS_TRANSPORT_H_


#include <iostream>
#include <vector>
#include <iterator>

#include "particle/particle.hpp"
#include "geometry/geometry.hpp"
#include "utils/rng.hpp"


enum BOUNDARY_CONDITION {
    REV = 0,  // periodic
    INJ = 1,  // re-injected
    OUT = 2   // out
};


class TransportMethod : public RNG
{

private: 

    // Systems

    std::shared_ptr<GeometrySystem> gs;
    std::shared_ptr<ParticleSystem> ps;

    // Files

    std::string fprefix;

    // Time parameters

    double t;
    double t0;
    double dt;
    double tmax;

    // Exit condition

    BOUNDARY_CONDITION bc;


public:

    std::vector<double> old_conc;

    TransportMethod(void);
    TransportMethod(int seed);
    TransportMethod(double _time_initial, double _time_delta, double _time_max, std::string _fprefix);
    TransportMethod(int seed, double _time_initial, double _time_delta, double _time_max, std::string _fprefix);

    // virtual void UpdateParticlePosition(int index_particle) = 0;

    // Setters

    void SetGeometrySystem(std::shared_ptr<GeometrySystem> _gs);
    void SetParticleSystem(std::shared_ptr<ParticleSystem> _ps);

    void SetFilePrefix(std::string _fprefix);

    void SetCurrentTime(double _time_current);
    void SetInitialTime(double _time_initial);
    void SetDeltaTime(double _time_delta);
    void SetMaximumTime(double _time_max);

    void SetBoundaryCondition(BOUNDARY_CONDITION _bc);

    // Getters

    int GetDimension(void);

    std::shared_ptr<GeometrySystem> GetGeometrySystem(void);
    std::shared_ptr<ParticleSystem> GetParticleSystem(void);

    std::string GetFilePrefix(void);

    double GetCurrentTime(void);
    double GetInitialTime(void);
    double GetDeltaTime(void);
    double GetMaximumTime(void);

    BOUNDARY_CONDITION GetBoundaryCondition(void);

    // Method

    std::vector<double> GetParticleConcentration(void);

    // I/O

    void WriteToFile(void);

};


//
// Constructors
//


inline TransportMethod::TransportMethod(void)
: RNG()
{}

inline TransportMethod::TransportMethod(int seed)
: RNG(seed)
{}

inline TransportMethod::TransportMethod(double _time_initial, double _time_delta, double _time_max, std::string _fprefix)
: RNG(), fprefix(_fprefix), t0(_time_initial), t(_time_initial), dt(_time_delta), tmax(_time_max) 
{}

inline TransportMethod::TransportMethod(int seed, double _time_initial, double _time_delta, double _time_max, std::string _fprefix)
: RNG(seed), fprefix(_fprefix), t0(_time_initial), t(_time_initial), dt(_time_delta), tmax(_time_max) 
{}

//
// Getters
//

inline int TransportMethod::GetDimension(void) 
{
    return this->GetGeometrySystem()->GetDimension();
}

inline std::shared_ptr<GeometrySystem> TransportMethod::GetGeometrySystem(void) 
{
    return gs;
}

inline std::shared_ptr<ParticleSystem> TransportMethod::GetParticleSystem(void) 
{
    return ps;
}

inline std::string TransportMethod::GetFilePrefix(void) 
{
    return fprefix;
}

inline double TransportMethod::GetCurrentTime(void) 
{
    return t;
}

inline double TransportMethod::GetInitialTime(void) 
{
    return t0;
}

inline double TransportMethod::GetDeltaTime(void) 
{
    return dt;
}

inline double TransportMethod::GetMaximumTime(void) 
{
    return tmax;
}

inline BOUNDARY_CONDITION TransportMethod::GetBoundaryCondition(void) 
{
    return bc;
}


//
// Setters
//

// Systems

inline void TransportMethod::SetFilePrefix(std::string _fprefix) 
{
    fprefix = _fprefix;
}

inline void TransportMethod::SetGeometrySystem(std::shared_ptr<GeometrySystem> _gs) 
{
    gs = std::move(_gs);
}

inline void TransportMethod::SetParticleSystem(std::shared_ptr<ParticleSystem> _ps) 
{
    ps = std::move(_ps);
}

// Time

inline void TransportMethod::SetCurrentTime(double _t) 
{
    t = _t;
}

inline void TransportMethod::SetBoundaryCondition(BOUNDARY_CONDITION _bc) 
{
    bc = _bc;
}

inline std::vector<double> TransportMethod::GetParticleConcentration(void) 
{
    for (int ip=0; ip<this->GetParticleSystem()->GetNumberOfParticles(); ++ip) 
    {
        int ind = sub2ind(
            this->GetParticleSystem()->GetGridSubIndex(ip), 
            this->GetGeometrySystem()->GetShape()
        );

        old_conc[ind]++;
    }
    
    return old_conc;
}

//
// I/O
//

inline void TransportMethod::WriteToFile(void) 
{
    ps->WriteParticleSummary((fprefix +                                             ".sum").c_str());
    ps->WriteParticleSystem((fprefix + "_" + std::to_string(ps->GetCurrentStep()) + ".vtp").c_str());
    // ToVTI((fprefix + "_C" + std::to_string(ps->GetCurrentStep())+".vti").c_str(),this->GetParticleConcentration(), this->GetGeometrySystem()->GetShape(), this->GetGeometrySystem()->GetSpacing(), this->GetGeometrySystem()->GetOrigin());
    // ps->WriteParticleProfile((fprefix + "_" + std::to_string(ps->GetCurrentStep()) + ".hist").c_str());
    // ps->WriteParticleXi((fprefix + "_" + std::to_string(ps->GetCurrentStep()) + ".hist").c_str(), v);
}

// Streams 

inline std::ostream &operator<<(std::ostream &os, const std::shared_ptr<TransportMethod> &tm) 
{
    os << message("Initial time",       tm->GetInitialTime());
    os << message("Delta time",         tm->GetDeltaTime()  );
    os << message("Max. physical time", tm->GetMaximumTime());

    return os;
}


#endif  // METHODS_TRANSPORT_H_
