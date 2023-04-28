
#ifndef METHODS_DEGRADE_HPP_
#define METHODS_DEGRADE_HPP_


#include <memory>

#include "utils/rng.hpp"
#include "methods/reactive.hpp"


class Degrade: public Reactive, RNG
{

private:

    // Constant

    double _k;
    double _dt;

public:

    Degrade(double k, double dt);

    void DegradeParticles(void);
    void UpdatePorosity(void);

};


inline Degrade::Degrade(double k, double dt)
: Reactive(), RNG(), _k(k), _dt(dt)
{}

inline void Degrade::DegradeParticles(void)
{
    for (int ip = 0; ip < this->GetParticleSystem()->GetNumberOfParticles(); ++ip) {
        if (this->GetParticleSystem()->GetState(ip) == 0) continue;

        double z = this->Sample();
        
        if (z < this->_k*this->_dt)
        {
            this->GetParticleSystem()->SetState(ip, 0);
        }
    }

}

inline void Degrade::UpdatePorosity(void)
{

    for (int ip = 0; ip < this->GetParticleSystem()->GetNumberOfParticles(); ++ip) {
        if (this->GetParticleSystem()->GetState(ip) == 1) continue;

        int index = this->GetGeometrySystem()->Sub2Ind(this->GetParticleSystem()->GetGridSubIndex(ip));
        double phi = 1.1*this->GetDarcyMethod()->GetPorosity(ip);

        if (phi >= 1)
        {
            phi = 1;
        }

        this->GetDarcyMethod()->SetPorosity(phi, index);
    }

}


#endif  // METHODS_DEGRADE_HPP_
