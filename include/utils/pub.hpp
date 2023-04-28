#ifndef PUB_HPP
#define PUB_HPP

#include <unordered_map>

#include "mpi.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkDataArrayCollection.h>
#include <vtkXMLPolyDataWriter.h>

#include "io/logger.hpp"
#include "utils/rng.hpp"
#include "geometry/geometry.hpp"
#include "particle/cell.hpp"
#include "particle/particle.hpp"
#include "methods/darcy.hpp"
#include "methods/psi.hpp"


enum INITIAL_CONDITION {
    POINT  = 0, 
    PLANE  = 1,
    VOLUME = 2,
    FLOW   = 3
};

enum BOUNDARY_CONDITION {
    PERIODIC = 0, 
    BT       = 1
};

std::unordered_map<std::string, INITIAL_CONDITION> map_initial_condition = {
    {"point",  INITIAL_CONDITION::POINT },
    {"plane",  INITIAL_CONDITION::PLANE },
    {"volume", INITIAL_CONDITION::VOLUME},
    {"flow",   INITIAL_CONDITION::FLOW  }
};

std::unordered_map<std::string, BOUNDARY_CONDITION> map_boundary_condition = {
    {"periodic",  BOUNDARY_CONDITION::PERIODIC},
    {"bt",        BOUNDARY_CONDITION::BT      }
};

class Pub : public RNG
{

public:

    int _seed = -1;

    int _step;
    double _t0, _t, _dt, _tmax;

    int num_species_fluid = 0;
    int num_species_solid = 0;
    int num_species_gas   = 0;
    int num_species_total = 0;

    int _bt = 0;
    Reaktoro::VectorXd _bt_species;

    std::vector<double> _x_range;
    std::vector<double> _y_range;
    std::vector<double> _z_range;

    std::vector<double> _mean;
    std::vector<double> _var;
    std::vector<double> _dl;
    std::vector<double> _var_old;

    int _nparts;
    int _nparts0;

    double _rho;
    double _molar; 
    double _rrate; 
    double _cmolar;

    std::string _fprefix;
    std::ofstream _fsummary;
    Logger *_log;

    MPI_Comm _comm;
    int _size;

    Cell _cells;
    std::shared_ptr<DarcyMethod> _darcy;

    ~Pub(void);
    Pub(void);
    Pub(MPI_Comm comm);

    void CreateParticleDistribution(int num_particles);
    void CreateParticleDistribution(int num_particles, double mols);
    void CreateParticleDistribution(int num_particles, Reaktoro::VectorXd amounts);
    void CreateParticleDistribution(int num_particles, Reaktoro::VectorXd amounts_0, Reaktoro::VectorXd amounts_1);
    
    void CreateParticleDistributionDisplacement(int num_particles, Reaktoro::VectorXd resident, Reaktoro::VectorXd injected);
    void CreateParticleUniformDistribution(void);

    void CreateParticlePlaneInjection(int num_particles, int ix);
    void CreateParticlePlaneInjection(int num_particles, int ix, double mols);
    void CreateParticlePlaneInjectionPhi0(int num_particles, int ix);
    void CreateParticlePlaneInjectionPhi0(int num_particles, int ix, double mols);
    void CreateParticlePlaneInjectionAB(int num_particles, int ix);
    // void CreateParticlePhiWeighted(int num_particles);
    void CreateParticlePointInjection(int num_particles, int ix, double mols);

    void CreateParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts);
    void CreateParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts_0, Reaktoro::VectorXd amounts_1);
    void CreateParticlePlaneInjectionPhi0(int num_particles, int ix, Reaktoro::VectorXd amounts);

    void AddParticlePlaneInjection(int num_particles, int ix);
    void AddParticlePlaneInjection(int num_particles, int ix, double mols);
    void AddParticlePlaneInjectionADE(int num_particles, int ix, double mols);
    void AddParticlePlaneInjectionPhi0(int num_particles, int ix, double mols);

    void AddParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts);
    void AddParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts_0, Reaktoro::VectorXd amounts_1);
    void AddParticlePlaneInjectionPhi0(int num_particles, int ix, Reaktoro::VectorXd amounts);

    double SampleTransitTime(int i, int j, int k, int dir);
    double SampleTransitTimeADE(int i, int j, int k, int dir);

    int SampleDirection(int i, int j, int k);

    void MoveNextLocation(std::list<Particle*>::iterator &p);

    void SetInitialConditions(void);
    void SetInitialConditionsADE(void);
    void SetTime(int step, double t0, double dt, double tmax);

    void IterateParticlesADE(void);
    void IterateParticlesADEDegrade(void);
    void IterateParticles(void);
    void IterateMixParticles(void);
    void IterateDegradeParticles(void);

    int GetNumberOfParticles(void);

    void CalculateParticleSummary(void);
    void CreateParticleSummaryFile(void);
    void WriteParticleSummary(void);

    void CreateSpeciesSummaryFile(ChemicalSystem system);
    void WriteSpeciesSummary(void);

    void WriteParticles(void);
    void ReadParticles(std::string fname);
    std::string OutputMessage(void);

};

Pub::~Pub(void)
{    
    MPI_Finalize();
}

Pub::Pub(void)
: RNG()
{}

Pub::Pub(MPI_Comm comm)
: _comm(comm), RNG()
{}

void Pub::CreateParticleDistribution(int num_particles)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> distx(1, n[0]-2);
    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int ix = distx(generator);
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::CreateParticleDistributionDisplacement(int num_particles, Reaktoro::VectorXd resident, Reaktoro::VectorXd injected) 
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int   > n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int   > c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> distx(2, n[0]-3);
    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int ix = distx(generator);
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        // if (ix == 2) // splitting the Y axis in 2 halves
        // {
        //     p->amounts = injected/static_cast<double>(num_particles);
        // }
        // else 
        // {
            p->amounts = resident/static_cast<double>(num_particles);
        // }

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::CreateParticleUniformDistribution(void) {
    std::vector<int>    n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);
    this->_cells = Cell(c, n);

    for (int k=0; k<n[2]; ++k) {
        for (int j=0; j<n[1]; ++j) {
            for (int i=0; i<n[0]; ++i) {
                int ind = sub2ind({i, j, k}, n);

                Particle* p = new Particle();
                p->_id = ind;
                p->_type = 0;
                p->_valid = true;

                p->SetX0(
                    o[0]+static_cast<double>(i)*s[0],
                    o[1]+static_cast<double>(j)*s[1],
                    o[2]+static_cast<double>(k)*s[2]
                );

                p->SetI0(i, j, k);
                p->SetT0(0.0, 0.0);

                p->SetX(p->_x0, p->_y0, p->_z0);
                p->SetI(p->_i0, p->_j0, p->_k0);
                
                this->_cells.push_back(p, i, j, k);
            }
        }
    }
}

void Pub::CreateParticleDistribution(int num_particles, Reaktoro::VectorXd amounts_0, Reaktoro::VectorXd amounts_1)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> distx(2, n[0]-3);
    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int ix = distx(generator);
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        if (iy <= n[1]/2) // splitting the Y axis in 2 halves
        {
            p->amounts = amounts_0/static_cast<double>(num_particles);
        }
        else 
        {
            p->amounts = amounts_1/static_cast<double>(num_particles);
        }

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::CreateParticleDistribution(int num_particles, Reaktoro::VectorXd amounts)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> distx(2, n[0]-3);
    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int ix = distx(generator);
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->amounts = amounts/static_cast<double>(num_particles);

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::CreateParticleDistribution(int num_particles, double mols)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> distx(2, n[0]-3);
    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int ix = distx(generator);
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;
        p->mols = mols/static_cast<double>(num_particles);

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}


void Pub::CreateParticlePlaneInjection(int num_particles, int ix)
{   
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;
        p->mols = 0;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::CreateParticlePointInjection(int num_particles, int ix, double mols)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = (n[1]-3)/2;
        int iz = (n[2]-3)/2;
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;
        p->mols = mols/static_cast<double>(num_particles);

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        this->_cells.push_back(p, ix, iy, iz); 
    }   
}

void Pub::CreateParticlePlaneInjection(int num_particles, int ix, double mols)
{   
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;
        p->mols = mols/static_cast<double>(num_particles);

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        this->_cells.push_back(p, ix, iy, iz); 
    }
}

void Pub::CreateParticlePlaneInjectionPhi0(int num_particles, int ix, double mols)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    this->_cells = Cell(c, n);

    xt::xarray<double> phi_plane = xt::view(this->_darcy->phi, ix, xt::all(), xt::all());
    double phi_max = xt::amax(phi_plane)[0];
    double phi_min = xt::amin(phi_plane)[0];

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        double phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);

        while (phi0 < 0.5)
        {
            iy = disty(generator);
            iz = distz(generator);
            ind = this->_darcy->Sub2Ind({ ix, iy, iz });

            phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);
        }

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;
        p->mols = mols/static_cast<double>(num_particles);

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::CreateParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts)
{   
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;
        p->amounts = amounts/static_cast<double>(num_particles);

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        this->_cells.push_back(p, ix, iy, iz); 
    }
}

void Pub::CreateParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts_0, Reaktoro::VectorXd amounts_1)
{   
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    this->_cells = Cell(c, n);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        if (iy <= n[1]/2) // splitting the Y axis in 2 halves
        {
            p->amounts = amounts_0/static_cast<double>(num_particles);
        }
        else 
        {
            p->amounts = amounts_1/static_cast<double>(num_particles);
        }

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        this->_cells.push_back(p, ix, iy, iz); 
    }
}

void Pub::CreateParticlePlaneInjectionPhi0(int num_particles, int ix, Reaktoro::VectorXd amounts)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::vector<int> c = std::vector<int>(3, 1);

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    this->_cells = Cell(c, n);

    xt::xarray<double> phi_plane = xt::view(this->_darcy->phi, ix, xt::all(), xt::all());
    double phi_max = xt::amax(phi_plane)[0];
    double phi_min = xt::amin(phi_plane)[0];

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);
        int ind = this->_darcy->Sub2Ind({ ix, iy, iz });

        double phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);

        while (phi0 < 0.5)
        {
            iy = disty(generator);
            iz = distz(generator);
            ind = this->_darcy->Sub2Ind({ ix, iy, iz });

            phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);
        }

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;
        p->amounts = amounts/static_cast<double>(num_particles);

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(0.0, 0.0);

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::AddParticlePlaneInjection(int num_particles, int ix)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(this->_t, 0.0);
        p->_step = this->_step;

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        p->_dir = this->SampleDirection(p->_i, p->_j, p->_k);
        p->_tt = this->SampleTransitTime(p->_i, p->_j, p->_k, p->_dir);
        p->_dtt = p->_tt;

        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::AddParticlePlaneInjectionADE(int num_particles, int ix, double mols)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(this->_t, 0.0);
        p->_step = this->_step;

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        p->_dir = this->SampleDirection(p->_i, p->_j, p->_k);
        p->_tt = this->SampleTransitTimeADE(p->_i, p->_j, p->_k, p->_dir);
        p->_dtt = p->_tt;

        p->mols = mols/static_cast<double>(num_particles);

        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::AddParticlePlaneInjection(int num_particles, int ix, double mols)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(this->_t, 0.0);
        p->_step = this->_step;

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        p->_dir = this->SampleDirection(p->_i, p->_j, p->_k);
        p->_tt = this->SampleTransitTime(p->_i, p->_j, p->_k, p->_dir);
        p->_dtt = p->_tt;

        p->mols = mols/static_cast<double>(num_particles);

        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::AddParticlePlaneInjectionPhi0(int num_particles, int ix, double mols)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    xt::xarray<double> phi_plane = xt::view(this->_darcy->phi, ix, xt::all(), xt::all());
    double phi_max = xt::amax(phi_plane)[0];
    double phi_min = xt::amin(phi_plane)[0];

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);

        double phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);

        while (phi0 < 0.5)
        {
            iy = disty(generator);
            iz = distz(generator);

            phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);
        }

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(this->_t, 0.0);
        p->_step = this->_step;

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        p->_dir = this->SampleDirection(p->_i, p->_j, p->_k);
        p->_tt = this->SampleTransitTime(p->_i, p->_j, p->_k, p->_dir);
        p->_dtt = p->_tt;

        p->mols = mols/static_cast<double>(num_particles);

        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::AddParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::uniform_int_distribution<int> disty(2, n[1]-3);
    std::uniform_int_distribution<int> distz(2, n[2]-3);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(this->_t, 0.0);
        p->_step = this->_step;

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        p->_dir = this->SampleDirection(p->_i, p->_j, p->_k);
        p->_tt = this->SampleTransitTime(p->_i, p->_j, p->_k, p->_dir);
        p->_dtt = p->_tt;

        p->amounts = amounts/static_cast<double>(num_particles);

        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::AddParticlePlaneInjection(int num_particles, int ix, Reaktoro::VectorXd amounts_0, Reaktoro::VectorXd amounts_1)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(this->_t, 0.0);
        p->_step = this->_step;

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        p->_dir = this->SampleDirection(p->_i, p->_j, p->_k);
        p->_tt = this->SampleTransitTime(p->_i, p->_j, p->_k, p->_dir);
        p->_dtt = p->_tt;

        if (iy <= n[1]/2) // splitting the Y axis in 2 halves
        {
            p->amounts = amounts_0/static_cast<double>(num_particles);
        }
        else 
        {
            p->amounts = amounts_1/static_cast<double>(num_particles);
        }
        
        this->_cells.push_back(p, ix, iy, iz);
    }
}

void Pub::AddParticlePlaneInjectionPhi0(int num_particles, int ix, Reaktoro::VectorXd amounts)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }

    std::vector<int> n = this->_darcy->_shape;
    std::vector<double> o = this->_darcy->_origin;
    std::vector<double> s = this->_darcy->_spacing;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    xt::xarray<double> phi_plane = xt::view(this->_darcy->phi, ix, xt::all(), xt::all());
    double phi_max = xt::amax(phi_plane)[0];
    double phi_min = xt::amin(phi_plane)[0];

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        int iy = disty(generator);
        int iz = distz(generator);

        double phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);

        while (phi0 < 0.5)
        {
            iy = disty(generator);
            iz = distz(generator);

            phi0 = (this->_darcy->phi(ix, iy, iz)-phi_min)/(phi_max-phi_min);
        }

        Particle* p = new Particle();

        p->_id = ip;

        p->_type = 0;
        p->_valid = true;

        p->SetX0(
            o[0]+static_cast<double>(ix)*s[0],
            o[1]+static_cast<double>(iy)*s[1],
            o[2]+static_cast<double>(iz)*s[2]
        );

        p->SetI0(ix, iy, iz);
        p->SetT0(this->_t, 0.0);
        p->_step = this->_step;

        p->SetX(p->_x0, p->_y0, p->_z0);
        p->SetI(p->_i0, p->_j0, p->_k0);

        p->_dir = this->SampleDirection(p->_i, p->_j, p->_k);
        p->_tt = this->SampleTransitTime(p->_i, p->_j, p->_k, p->_dir);
        p->_dtt = p->_tt;

        p->amounts = amounts/static_cast<double>(num_particles);

        this->_cells.push_back(p, ix, iy, iz);
    }
}

double Pub::SampleTransitTime(int i, int j, int k, int dir)
{
    return this->_darcy->psi(i, j, k, dir).SampleTime(this->Sample());
}

double Pub::SampleTransitTimeADE(int i, int j, int k, int dir)
{
    return this->_darcy->psi_ade(i, j, k).SampleTime(this->Sample(), dir);
}

int Pub::SampleDirection(int i, int j, int k)
{
    xt::xarray<double> temp_prob = 
        xt::view(this->_darcy->direction_prob_cumm, i, j, k, xt::range(0, 7));

    int ind = 0;
    double z = this->Sample();
    
    for (double v : temp_prob)
    {
        if (v > z) break;
        ++ind;
    }

    if (ind >= 6)
    {
        return 5;
    }
    else if (ind <= 0)
    {
        return 0;
    }
    else
    {
        return ind;
    }
}

void Pub::SetInitialConditions(void)
{
    this->_mean    = std::vector<double>(3, 0);
    this->_var     = std::vector<double>(3, 0);
    this->_var_old = std::vector<double>(3, 0);
    this->_dl      = std::vector<double>(3, 0);

    for (std::list<Particle*> plist : this->_cells._parts)
    {
        std::list<Particle*>::iterator p = plist.begin();
    
        while (p != plist.end())
        {
            (*p)->_dir = this->SampleDirection((*p)->_i, (*p)->_j, (*p)->_k);
            (*p)->_tt  = this->SampleTransitTime((*p)->_i, (*p)->_j, (*p)->_k, (*p)->_dir);
            (*p)->_dtt = (*p)->_tt;

            (*p)->v = this->_darcy->CalculateDistance((*p)->GetSourceI(), (*p)->GetTargetI())/(*p)->_tt;

            p++;
        }
    }
}

void Pub::SetInitialConditionsADE(void)
{
    this->_mean    = std::vector<double>(3, 0);
    this->_var     = std::vector<double>(3, 0);
    this->_var_old = std::vector<double>(3, 0);
    this->_dl      = std::vector<double>(3, 0);

    for (std::list<Particle*> plist : this->_cells._parts) {
        std::list<Particle*>::iterator p = plist.begin();
    
        while (p != plist.end()) {
            (*p)->_dir = this->SampleDirection((*p)->_i, (*p)->_j, (*p)->_k);
            (*p)->_tt = this->SampleTransitTimeADE((*p)->_i, (*p)->_j, (*p)->_k, (*p)->_dir);
            (*p)->_dtt = (*p)->_tt;

            (*p)->v = this->_darcy->CalculateDistance((*p)->GetSourceI(), (*p)->GetTargetI())/(*p)->_tt;

            p++;
        }
    }
}

void Pub::SetTime(int step, double t0, double dt, double tmax)
{
    this->_step = step;
    this->_t0   = t0;
    this->_t    = t0;
    this->_dt   = dt;
    this->_tmax = tmax;

    for (std::list<Particle*> plist : this->_cells._parts)
    {
        std::list<Particle*>::iterator p = plist.begin();

        while (p != plist.end())
        {
            (*p)->SetT0(t0, 0.0);
            (*p)->_step = step;
            ++p;
        }
    }
}

void Pub::IterateParticlesADE(void)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }
    
    std::vector<int> n = this->_darcy->_shape;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    this->CreateParticleSummaryFile();

    double dt = this->_dt;
    double tmax = this->_tmax;

    while (this->_t < tmax)
    {
        for (std::list<Particle*> &plist : this->_cells._parts)
        {
            std::list<Particle*>::iterator p = plist.begin();

            while (p != plist.end() && (*p)->_valid && (*p)->_t < this->_t)
            {
                if ((*p)->_dtt > dt)
                {
                    (*p)->_t += (*p)->_dtt;
                    (*p)->_dtt = 0.0f;
                    (*p)->_move_list = false;
                }
                else
                {
                    xt::xarray<int> dir = xt::view(direction_vector(), (*p)->_dir, xt::all());
                    xt::xarray<int> ind { (*p)->_i, (*p)->_j, (*p)->_k };

                    ind += dir;

                    if (this->_darcy->IsBoundary(ind(0), ind(1), ind(2)))
                    {
                        if (this->_darcy->IsXp(ind))
                        {
                            /*
                            if (REVX == true)
                            {
                                ++(this->_bt);
                                // (*p)->_valid = false;
                                ind(0) = 1;
                                ind(1) = disty(generator);
                                ind(2) = distz(generator);
                            }
                            */
                            ++(this->_bt);
                            (*p)->_valid = false;
                        }
                        else
                        {
                            ind = this->_darcy->GetBouncedGridSubindex(ind);
                        }
                    }

                    if ((*p)->_valid)
                    {
                        std::vector<double> x = this->_darcy->GetLocation((*p)->_i, (*p)->_j, (*p)->_k);

                        (*p)->SetX(x.at(0), x.at(1), x.at(2));
                        (*p)->SetI(ind(0), ind(1), ind(2));

                        (*p)->_dir = this->SampleDirection(ind(0), ind(1), ind(2));
                        (*p)->_tt = this->SampleTransitTimeADE(ind(0), ind(1), ind(2), (*p)->_dir);
                        (*p)->_dtt = (*p)->_dtt + (*p)->_tt;

                        (*p)->_move_list = true;
                    }
                }

                if (!(*p)->_valid)
                {
                    p = plist.erase(p);
                }
                else if ((*p)->_move_list)
                {
                    (*p)->_move_list = false;

                    std::list<Particle*> &plist_target = this->_cells._parts((*p)->_i,  (*p)->_j,  (*p)->_k );
                    (*p)->SetI0((*p)->_i, (*p)->_j, (*p)->_k);
                    
                    plist_target.splice(plist_target.begin(), plist, p++);
                }
                else
                {
                    ++p;
                }
            }
        }

        this->_step += 1;
        this->_t += dt;

        this->CalculateParticleSummary();
        this->WriteParticleSummary();
        this->WriteParticles();

        *(this->_log) << this->OutputMessage() << std::endl;
    }
}

void Pub::IterateParticlesADEDegrade(void)
{
    std::random_device device;
    std::mt19937 generator(device());

    if (this->_seed =! -1)
    {
        generator.seed(this->_seed);
    }
    
    std::vector<int> n = this->_darcy->_shape;

    std::uniform_int_distribution<int> disty(1, n[1]-2);
    std::uniform_int_distribution<int> distz(1, n[2]-2);

    this->CreateParticleSummaryFile();

    double dt = this->_dt;
    double tmax = this->_tmax;

    while (this->_t < tmax)
    {
        for (std::list<Particle*> &plist : this->_cells._parts)
        {
            std::list<Particle*>::iterator p = plist.begin();

            while (p != plist.end() && (*p)->_valid && (*p)->_t < this->_t)
            {
                if ((*p)->_dtt > dt)
                {
                    (*p)->_t += (*p)->_dtt;
                    (*p)->_dtt = 0.0f;
                    (*p)->_move_list = false;
                }
                else
                {
                    xt::xarray<int> dir = xt::view(direction_vector(), (*p)->_dir, xt::all());
                    xt::xarray<int> ind { (*p)->_i, (*p)->_j, (*p)->_k };

                    ind += dir;

                    if (this->_darcy->IsBoundary(ind(0), ind(1), ind(2)))
                    {
                        if (this->_darcy->IsXp(ind))
                        {
                            ++(this->_bt);
                            (*p)->_valid = false;
                        }
                        else
                        {
                            ind = this->_darcy->GetBouncedGridSubindex(ind);
                        }
                    }

                    if ((*p)->_valid)
                    {
                        std::vector<double> x = this->_darcy->GetLocation((*p)->_i, (*p)->_j, (*p)->_k);

                        (*p)->SetX(x.at(0), x.at(1), x.at(2));
                        (*p)->SetI(ind(0), ind(1), ind(2));

                        (*p)->_dir = this->SampleDirection(ind(0), ind(1), ind(2));
                        (*p)->_tt = this->SampleTransitTimeADE(ind(0), ind(1), ind(2), (*p)->_dir);
                        (*p)->_dtt = (*p)->_dtt + (*p)->_tt;

                        (*p)->_move_list = true;
                    }
                }

                if (!(*p)->_valid)
                {
                    p = plist.erase(p);
                }
                else if ((*p)->_move_list)
                {
                    (*p)->_move_list = false;

                    std::list<Particle*> &plist_target = this->_cells._parts((*p)->_i,  (*p)->_j,  (*p)->_k );
                    (*p)->SetI0((*p)->_i, (*p)->_j, (*p)->_k);
                    
                    plist_target.splice(plist_target.begin(), plist, p++);
                }
                else
                {
                    ++p;
                }
            }
        }

        for (std::list<Particle*> &plist : this->_cells._parts)
        {
            std::list<Particle*>::iterator p = plist.begin();

            int N = 0;
        
            while (p != plist.end())
            {
                ++N;
                ++p;
            }

            if (N > 0)
            {
                std::list<Particle*>::iterator p = plist.begin();

                int i = (*p)->_i;
                int j = (*p)->_j;
                int k = (*p)->_k;

                double f = 1.1;
                double kr = this->_darcy->phi(i, j, k);
                if (kr > 1)
                {
                    kr = 0;
                }
                else 
                {
                    kr = 1-kr;
                }

                int n = 0;

                if (kr > 0.0001)
                {
                    double r = this->Sample();
                    double tau = std::log(1/r)/(N*kr);
                    n = static_cast<int>(dt/tau);
                }

                while (n > 0 && p != plist.end())
                {
                    this->_darcy->phi(i, j, k) = f*this->_darcy->phi(i, j, k);

                    if (this->_darcy->phi(i, j, k) > 1.0)
                    {
                        this->_darcy->phi(i, j, k) = 1.0;
                    }

                    --n;
                    p = plist.erase(p);
                    // ++p;
                }

            }
        }

        this->_step += 1;
        this->_t += dt;

        if (std::remainder(this->_step, 5.0) == 0)
        {
            *(this->_log) << message("Writing Pressure in step", this->_step) << std::endl;
            this->_darcy->CalculateCorrelatedPermeability(0.0001, 10.0);
            this->_darcy->SolvePressure(this->_comm, 0, 0);
            this->_darcy->CalculateVelocity();
            this->_darcy->CalculateProbabilityADE();

            this->_darcy->WritePressureField((this->_fprefix+"_pressure"+std::to_string(this->_step)+".vti").c_str());
            xt::dump_npy((this->_fprefix+"_pressure"+std::to_string(this->_step)+".npy").c_str(), this->_darcy->pressure);

            this->_darcy->WritePermeabilityField((this->_fprefix+"_perm"+std::to_string(this->_step)+".vti").c_str());
            xt::dump_npy((this->_fprefix+"_perm"+std::to_string(this->_step)+".npy").c_str(), this->_darcy->perm);
            
            this->_darcy->WritePorosityField((this->_fprefix+"_poro"+std::to_string(this->_step)+".vti").c_str());
            xt::dump_npy((this->_fprefix+"_poro"+std::to_string(this->_step)+".npy").c_str(), this->_darcy->phi);
        }

        this->CalculateParticleSummary();
        this->WriteParticleSummary();
        this->WriteParticles();

        *(this->_log) << this->OutputMessage() << std::endl;
        this->AddParticlePlaneInjection(this->_nparts0, 2);
    }
}

void Pub::IterateParticles(void)
{
    this->CreateParticleSummaryFile();

    double dt = this->_dt;
    double tmax = this->_tmax;

    while (this->_t < tmax)
    {
        for (std::list<Particle*> &plist : this->_cells._parts)
        {
            std::list<Particle*>::iterator p = plist.begin();

            while (p != plist.end() && (*p)->_t < this->_t)
            {
                if ((*p)->_dtt > dt)
                {
                    (*p)->_t += (*p)->_dtt;
                    (*p)->_dtt = 0.0f;
                    (*p)->_move_list = false;
                }
                else
                {
                    // this->MoveNextLocation(p);

                    xt::xarray<int> dir = xt::view(direction_vector(), (*p)->_dir, xt::all());
                    xt::xarray<double> x { (*p)->_x, (*p)->_y, (*p)->_z };
                    xt::xarray<int> ind { (*p)->_i, (*p)->_j, (*p)->_k };

                    x += xt::adapt(this->_darcy->_spacing, {3})*dir;
                    ind += dir;

                    (*p)->SetX(x(0), x(1), x(2));
                    (*p)->SetI(ind(0), ind(1), ind(2));

                    (*p)->_dir = this->SampleDirection(ind(0), ind(1), ind(2));
                    (*p)->_tt = this->SampleTransitTime(ind(0), ind(1), ind(2), (*p)->_dir);
                    (*p)->_dtt = (*p)->_dtt + (*p)->_tt;

                    if (this->_darcy->IsBoundary(ind(0), ind(1), ind(2)))
                    {
                        (*p)->_valid = false;
                    }

                    (*p)->_move_list = true;    
                }

                if (!(*p)->_valid)
                {
                    p = plist.erase(p);
                }
                else if ((*p)->_move_list)
                {
                    (*p)->_move_list = false;

                    std::list<Particle*> &plist_target = this->_cells._parts((*p)->_i,  (*p)->_j,  (*p)->_k );
                    (*p)->SetI0((*p)->_i, (*p)->_j, (*p)->_k);
                    
                    plist_target.splice(plist_target.begin(), plist, p++);
                }
                else
                {
                    ++p;
                }
            }
        }

        this->_step += 1;
        this->_t += dt;

        this->CalculateParticleSummary();
        this->WriteParticleSummary();
        this->WriteParticles();

        *(this->_log) << this->OutputMessage() << std::endl;
    }
}
/*
void Pub::IterateMixParticles(void)
{
    this->CreateParticleSummaryFile();

    double dt = this->_data.GetValue<double>("dt");
    double tmax = this->_data.GetValue<double>("tmax");

    while (this->_data.GetValue<double>("t") < tmax)
    {
        // CTRW

        for (Particle* head : this->_cells._parts)
        {
            Particle* p = head;

            while (p != nullptr)
            {
                while (p->_data.GetValue<double>("t") < this->_data.GetValue<double>("t")) 
                {
                    if (p->_data.GetValue<double>("dtt") > dt)
                    {
                        p->_data.IncrementValue<double>("t", p->_data.GetValue<double>("dtt"));
                        p->_data.SetValue("dtt", 0.0f);
                        // p->_data.IncrementValue<double>("dtt", this->_data.GetValue<double>("t")+dt-p->GetElapsedT());
                    }
                    else
                    {
                        // p->_data.IncrementValue<double>("t", p->GetTargetT()-p->GetElapsedT());
                        // p->_data.SetValue("dtt", 0.0f);

                        this->MoveNextLocation(p);
                        p->_data.SetValue("dir", this->SampleDirection(p));

                        p->_data.SetValue(
                            "tt",
                            this->SampleTransitTime(
                                p->_i,
                                p->_j,
                                p->_k,
                                p->_dir
                            )
                        );

                        p->_data.IncrementValue<double>("dtt", p->_data.GetValue<double>("tt"));
                    }

                    p->_data.IncrementValue<int>("step", 1);
                    // TODO: correct particle for half-travel?
                }

                p = p->next;

                // if (!p->prev->_data.GetValue<bool>("valid"))
                // {
                //     this->_cells.remove(head, p);
                // }
            }
        }

        // Mixing

        for (Particle* head : this->_cells._parts)
        {
            Particle* p = head;

            int countA = 0;
            int countB = 0;

            while (p != nullptr)
            {
                if (p->_data.GetValue<int>("type") == 0)
                {
                    ++countA;
                }
                else
                {
                    ++countB;
                }

                p = p->next;
            }

            if (countA > 0 && countB > 0)
                std::cout << " A: " << countA << " B: " << countB << std::endl;

            double k = 0.1f;
            double A = static_cast<double>(countA)/static_cast<double>(countA+countB);
            double B = static_cast<double>(countB)/static_cast<double>(countA+countB);

            if (countA && countB && this->Sample() > k*A*B*dt)
            {
                bool removed = false;
                bool converted = false;

                p = head;

                while (p != nullptr)
                {
                    if (!removed && p->_data.GetValue<int>("type") == 0)
                    {
                        this->_cells.remove(p);
                        removed = true;
                    }

                    if (!converted && p->_data.GetValue<int>("type") == 1)
                    {
                        p->_data.SetValue("type", 2);
                        converted = true;
                    }

                    p = p->next;
                }   
            }
        }

        this->_data.IncrementValue<int>("step", 1);
        this->_data.IncrementValue<double>("t", dt);

        this->CalculateParticleSummary();
        this->WriteParticleSummary();
        this->WriteParticles();

        *(this->_log) << this->OutputMessage() << std::endl;
    }
}
*/
/*
void Pub::IterateDegradeParticles(void)
{
    this->CreateParticleSummaryFile();

    double dt = this->_dt;
    double tmax = this->_tmax;

    while (this->_t < tmax)
    {
        // CTRW
        for (std::list<Particle*> plist : this->_cells._parts)
        {
            std::list<Particle*>::iterator p = plist.begin();

            while (p != plist.end())
            {
                while (p->_t < this->_t && p->_valid) 
                {
                    if (p->_dtt > dt)
                    {
                        p->_t = p->_t + p->_dtt;
                        p->_dtt = 0.0f;
                    }
                    else
                    {
                        std::cout << " == Move == " << std::endl;
                        this->MoveNextLocation(p);
                        std::cout << " == p_ij == " << std::endl;
                        p->_dir = this->SampleDirection(p);
                        std::cout << " == psi == " << std::endl;
                        p->_tt = this->SampleTransitTime(p->_i, p->_j, p->_k, p->_dir);

                        p->_dtt = p->_dtt + p->_tt;
                    }

                    p->_step = p->_step + 1;
                    // TODO: correct particle for half-travel?
                }

                ++p;
            }
        }

        // Mixing
        for (std::list<Particle*> plist : this->_cells._parts)
        {
            std::list<Particle*>::iterator p = plist.begin();

            int N = 0;
            
            while (p != plist.end())
            {
                ++N;
                ++p;
            }

            if (N > 0) 
            {
                p = head;

                int i = p->_i;
                int j = p->_j;
                int k = p->_k;

                double f = 0.5f;
                double kr = 1*this->_darcy->propensity(i, j, k);
                int n = 0;

                if (kr > 0.0001)
                {
                    double r = this->Sample();
                    double tau = std::log(1/r)/(N*kr);
                    n = static_cast<int>(dt/tau);
                }

                while (n > 0 && p != nullptr)
                {
                    Particle* temp = this->_cells.remove(p);
                    
                    this->_darcy->phi(i, j, k) = f*this->_darcy->phi(i, j, k);

                    if (this->_darcy->phi(i, j, k) < 0.01)
                    {
                        this->_darcy->phi(i, j, k) = 0.0;
                    }
                    
                    --n;
                    ++p;
                }
            }
        }

        this->_step = 1;
        this->_t = dt;

        this->CalculateParticleSummary();
        this->WriteParticleSummary();
        
        if (std::remainder(this->_step, 5.0) == 0)
        {
            std::cout << "saving pressure in " << this->_step << std::endl;
            this->_darcy->SolvePressure(this->_comm, 0, 0);
            this->_darcy->WritePorosityField(
                (this->_fprefix+"_poro"+std::to_string(this->_step)+".vti").c_str()
            );
        }   

        this->WriteParticles();

        *(this->_log) << this->OutputMessage() << std::endl;
        this->AddParticlePlaneInjection(this->_nparts0, 1);
    }
}
*/
int Pub::GetNumberOfParticles(void)
{
    int n = 0;

    for (std::list<Particle*> plist : this->_cells._parts)
    {
        n += plist.size();
    }

    this->_nparts = n;
    
    return n;
}

void Pub::CalculateParticleSummary(void)
{
    int n = this->GetNumberOfParticles();

    this->_x_range = std::vector<double>({
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::lowest()
    });

    this->_y_range = std::vector<double>({
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::lowest()
    });

    this->_z_range = std::vector<double>({
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::lowest()
    });

    for (std::list<Particle*> &plist : this->_cells._parts) {
        std::list<Particle*>::iterator p = plist.begin();

        while (p != plist.end()) {
            if ((*p)->_x < this->_x_range[0]) this->_x_range[0] = (*p)->_x;
            if ((*p)->_x > this->_x_range[1]) this->_x_range[1] = (*p)->_x;

            if ((*p)->_y < this->_y_range[0]) this->_y_range[0] = (*p)->_y;
            if ((*p)->_y > this->_y_range[1]) this->_y_range[1] = (*p)->_y;

            if ((*p)->_z < this->_z_range[0]) this->_z_range[0] = (*p)->_z;
            if ((*p)->_z > this->_z_range[1]) this->_z_range[1] = (*p)->_z;

            this->_mean[0] += (*p)->_x;
            this->_mean[1] += (*p)->_y;
            this->_mean[2] += (*p)->_z;

            ++p;
        }
    }

    for (int i=0; i<3; ++i) this->_mean[i] /= static_cast<double>(n);

    for (std::list<Particle*> &plist : this->_cells._parts) {
        std::list<Particle*>::iterator p = plist.begin();

        while (p != plist.end()) {
            this->_var[0] += std::pow((*p)->_x-this->_mean[0], 2);
            this->_var[1] += std::pow((*p)->_y-this->_mean[1], 2);
            this->_var[2] += std::pow((*p)->_z-this->_mean[2], 2);
            
            ++p;
        }
    }

    for (int i=0; i<3; ++i) this->_var[i] /= static_cast<double>(n);
    for (int i=0; i<3; ++i) this->_dl[i] = (this->_var[i] - this->_var_old[i]) / this->_dt;

    this->_var_old = this->_var;

}

void Pub::CreateParticleSummaryFile(void)
{
    this->_fsummary.open(
        this->_fprefix + ".sum",
        std::ofstream::out
    );

    this->_fsummary << std::right << "step,t,";
    this->_fsummary << "x_min,y_min,z_min,";
    this->_fsummary << "x_max,y_max,z_max,";
    this->_fsummary << "x_mean,y_mean,z_mean,";
    this->_fsummary << "x_var,y_var,z_var,";
    this->_fsummary << "x_dl,y_dl,z_dl,";
    this->_fsummary << "poutflow" << std::endl;

    this->_fsummary.close();
}

void Pub::WriteParticleSummary(void)
{
    this->_fsummary.open(
        this->_fprefix + ".sum", 
        std::ofstream::out | std::ofstream::app
    );
    
    this->_fsummary << std::right << this->_step << ",";
    this->_fsummary << this->_t << ",";

    this->_fsummary << this->_x_range[0] << ",";
    this->_fsummary << this->_y_range[0] << ",";
    this->_fsummary << this->_z_range[0] << ",";
    this->_fsummary << this->_x_range[1] << ",";
    this->_fsummary << this->_y_range[1] << ",";
    this->_fsummary << this->_z_range[1] << ",";

    this->_fsummary << this->_mean[0] << ",";
    this->_fsummary << this->_mean[1] << ",";
    this->_fsummary << this->_mean[2] << ",";
    
    this->_fsummary << this->_var[0] << ",";
    this->_fsummary << this->_var[1] << ",";
    this->_fsummary << this->_var[2] << ",";
    
    this->_fsummary << this->_dl[0] << ",";
    this->_fsummary << this->_dl[1] << ",";
    this->_fsummary << this->_dl[2] << ",";

    this->_fsummary << this->_bt << std::endl;

    this->_fsummary.close();
}

void Pub::CreateSpeciesSummaryFile(ChemicalSystem system)
{
    this->_fsummary.open(
        this->_fprefix + ".species",
        std::ofstream::out
    );

    this->_fsummary << std::right << "step,t,";

    for (int i=0; i<this->num_species_fluid-1; ++i) {
        this->_fsummary << system.species(i).name() << ",";
    }

    this->_fsummary << system.species(this->num_species_fluid-1).name() << std::endl;
    this->_fsummary.close();
}

void Pub::WriteSpeciesSummary(void)
{
    this->_fsummary.open(
        this->_fprefix + ".species",
        std::ofstream::out | std::ofstream::app
    );

    this->_fsummary << std::right << this->_step << ",";
    this->_fsummary << this->_t << ",";

    for (int i=0; i<this->num_species_fluid-1; ++i) {
        this->_fsummary << this->_bt_species(i) << ",";
    }

    this->_fsummary << this->_bt_species(this->num_species_fluid-1) << std::endl;
    this->_fsummary.close();
}

void Pub::WriteParticles(void) {
    int n = this->GetNumberOfParticles();

    // Point data collection

    vtkSmartPointer<vtkDataArrayCollection> point_collection =
        vtkSmartPointer<vtkDataArrayCollection>::New();

    // Field data collection

    vtkSmartPointer<vtkDataArrayCollection> field_collection = 
        vtkSmartPointer<vtkDataArrayCollection>::New();

    // Point data arrays

    vtkSmartPointer<vtkIntArray> id_array = 
        vtkSmartPointer<vtkIntArray>::New();
    id_array->SetName("id");
    id_array->SetNumberOfComponents(1);
    id_array->Allocate(n);
    point_collection->AddItem(id_array);

    vtkSmartPointer<vtkIntArray> move_array = 
        vtkSmartPointer<vtkIntArray>::New();
    move_array->SetName("move");
    move_array->SetNumberOfComponents(1);
    move_array->Allocate(n);
    point_collection->AddItem(move_array);

    vtkSmartPointer<vtkDoubleArray> x0_array =
        vtkSmartPointer<vtkDoubleArray>::New();
    x0_array->SetName("x0");
    x0_array->SetNumberOfComponents(3);
    x0_array->Allocate(n);
    point_collection->AddItem(x0_array);

    vtkSmartPointer<vtkDoubleArray> xg_array =
        vtkSmartPointer<vtkDoubleArray>::New();
    xg_array->SetName("xg");
    xg_array->SetNumberOfComponents(3);
    xg_array->Allocate(n);
    point_collection->AddItem(xg_array);

    vtkSmartPointer<vtkIntArray> index0_array =
        vtkSmartPointer<vtkIntArray>::New();
    index0_array->SetName("index0_array");
    index0_array->SetNumberOfComponents(3);
    index0_array->Allocate(n);
    point_collection->AddItem(index0_array);

    vtkSmartPointer<vtkIntArray> source_array =
        vtkSmartPointer<vtkIntArray>::New();
    source_array->SetName("source");
    source_array->SetNumberOfComponents(3);
    source_array->Allocate(n);
    point_collection->AddItem(source_array);

    vtkSmartPointer<vtkIntArray> target_array =
        vtkSmartPointer<vtkIntArray>::New();
    target_array->SetName("target");
    target_array->SetNumberOfComponents(3);
    target_array->Allocate(n);
    point_collection->AddItem(target_array);

    vtkSmartPointer<vtkIntArray> dir_array = 
        vtkSmartPointer<vtkIntArray>::New();
    dir_array->SetName("dir");
    dir_array->SetNumberOfComponents(1);
    dir_array->Allocate(n);
    point_collection->AddItem(dir_array);

    vtkSmartPointer<vtkDoubleArray> t0_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    t0_array->SetName("t0");
    t0_array->SetNumberOfComponents(1);
    t0_array->Allocate(n);
    point_collection->AddItem(t0_array);

    vtkSmartPointer<vtkDoubleArray> t_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    t_array->SetName("t");
    t_array->SetNumberOfComponents(1);
    t_array->Allocate(n);
    point_collection->AddItem(t_array);

    vtkSmartPointer<vtkDoubleArray> tt_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    tt_array->SetName("tt");
    tt_array->SetNumberOfComponents(1);
    tt_array->Allocate(n);
    point_collection->AddItem(tt_array);

    vtkSmartPointer<vtkDoubleArray> dtt_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    dtt_array->SetName("dtt");
    dtt_array->SetNumberOfComponents(1);
    dtt_array->Allocate(n);
    point_collection->AddItem(dtt_array);

    vtkSmartPointer<vtkDoubleArray> v_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    v_array->SetName("v");
    v_array->SetNumberOfComponents(1);
    v_array->Allocate(n);
    point_collection->AddItem(v_array);

    vtkSmartPointer<vtkDoubleArray> mol_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    mol_array->SetName("mol");
    mol_array->SetNumberOfComponents(1);
    mol_array->Allocate(n);
    point_collection->AddItem(mol_array);

    vtkSmartPointer<vtkIntArray> type_array = 
        vtkSmartPointer<vtkIntArray>::New();
    type_array->SetName("type");
    type_array->SetNumberOfComponents(1);
    type_array->Allocate(n);
    point_collection->AddItem(type_array);

    vtkSmartPointer<vtkIntArray> valid_array = 
        vtkSmartPointer<vtkIntArray>::New();
    valid_array->SetName("valid");
    valid_array->SetNumberOfComponents(1);
    valid_array->Allocate(n);
    point_collection->AddItem(valid_array);

    vtkSmartPointer<vtkDoubleArray> species_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    
    if (this->num_species_fluid > 0) {
        species_array->SetName("species");
        species_array->SetNumberOfComponents(this->num_species_fluid);
        species_array->Allocate(n);
        point_collection->AddItem(species_array);
    }
    
    // Field data arrays

    vtkSmartPointer<vtkIntArray> step_array = 
        vtkSmartPointer<vtkIntArray>::New();
    step_array->SetName("step");
    step_array->SetNumberOfComponents(1);
    field_collection->AddItem(step_array);

    vtkSmartPointer<vtkDoubleArray> tg_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    tg_array->SetName("t");
    tg_array->SetNumberOfComponents(1);
    field_collection->AddItem(tg_array);

    vtkSmartPointer<vtkDoubleArray> x_mean_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    x_mean_array->SetName("x_mean");
    x_mean_array->SetNumberOfComponents(3);
    field_collection->AddItem(x_mean_array);

    vtkSmartPointer<vtkDoubleArray> x_var_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    x_var_array->SetName("x_var");
    x_var_array->SetNumberOfComponents(3);
    field_collection->AddItem(x_var_array);

    vtkSmartPointer<vtkDoubleArray> x_dl_array = 
        vtkSmartPointer<vtkDoubleArray>::New();
    x_dl_array->SetName("x_dl");
    x_dl_array->SetNumberOfComponents(3);
    field_collection->AddItem(x_dl_array);

    vtkSmartPointer<vtkIntArray> poutflow_array = 
        vtkSmartPointer<vtkIntArray>::New();
    poutflow_array->SetName("bt");
    poutflow_array->SetNumberOfComponents(1);
    field_collection->AddItem(poutflow_array);

    // Points

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    points->Allocate(n);

    // Cells

    vtkSmartPointer<vtkCellArray> verts =
        vtkSmartPointer<vtkCellArray>::New();

    // Set data array values from nodes
    int count = 0;

    for (std::list<Particle*> plist : this->_cells._parts)
    {
        std::list<Particle*>::iterator p = plist.begin();

        while (p != plist.end())
        {
            vtkIdType vtkid[1] = { count };  // Confused why you can't use ints or
                                             // unsigned ints in some places?
                                             // http://www.vtk.org/Wiki/VTK/Tutorials/VtkIdType

            std::vector<double> x = { (*p)->_x, (*p)->_y, (*p)->_z };
            xt::xarray<int> source = (*p)->GetSourceI();
            xt::xarray<int> target = (*p)->GetTargetI();

            points->InsertNextPoint(x.data());
            verts->InsertNextCell(1, vtkid);

            xg_array->InsertNextTuple3((*p)->xg[0], (*p)->xg[1], (*p)->xg[2]);

            index0_array->InsertNextTuple3((*p)->_i0, (*p)->_j0, (*p)->_k0);
            source_array->InsertNextTuple3(source(0), source(1), source(2));
            target_array->InsertNextTuple3(target(0), target(1), target(2));
            dir_array->InsertNextValue((*p)->_dir);

            id_array->InsertNextValue((*p)->_id);
            move_array->InsertNextValue((*p)->_move_list);
            type_array->InsertNextValue((*p)->_type);
            valid_array->InsertNextValue((*p)->_valid);

            x0_array->InsertNextTuple3((*p)->_x0, (*p)->_y0, (*p)->_z0);
            
            t0_array->InsertNextValue((*p)->_t0);
            t_array->InsertNextValue((*p)->_t);
            tt_array->InsertNextValue((*p)->_tt);
            dtt_array->InsertNextValue((*p)->_tt);

            v_array->InsertNextValue((*p)->v);

            mol_array->InsertNextValue((*p)->mols);

            if (this->num_species_fluid > 0) {
                species_array->InsertNextTuple((*p)->amounts.data());
            }

            ++count;
            ++p;
        }
    }

    step_array->InsertNextValue(this->_step);
    tg_array->InsertNextValue(this->_t);
    x_mean_array->InsertNextTuple(this->_mean.data());
    x_var_array->InsertNextTuple(this->_var.data());
    x_dl_array->InsertNextTuple(this->_dl.data());
    poutflow_array->InsertNextValue(this->_bt);

    // PolyData dataset

    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();

    polydata->SetPoints(points);
    polydata->SetVerts (verts);

    for (int i = 0; i < point_collection->GetNumberOfItems(); ++i)
        polydata->GetPointData()->AddArray(point_collection->GetItem(i));

    for (int i = 0; i < field_collection->GetNumberOfItems(); ++i)
        polydata->GetFieldData()->AddArray(field_collection->GetItem(i));   

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(
        (
            this->_fprefix + "_" + std::to_string(this->_step) + ".vtp"
        ).c_str()
    );
    writer->SetInputData(polydata);
    writer->Write();

}

void Pub::ReadParticles(std::string fname) {
    
}

std::string Pub::OutputMessage(void)
{
    std::stringstream ss;

    ss << message("Number of particles", this->_nparts);
    ss << message("Current step", this->_step);
    ss << message("Elapsed time [s]", this->_t);

    ss << message("X Range [m]", this->_x_range);
    ss << message("Y Range [m]", this->_y_range);
    ss << message("Z Range [m]", this->_z_range);

    ss << message("Mean position [m]", this->_mean);
    ss << message("Var position", this->_var);
    ss << message("Dl", this->_dl);

    return ss.str();
}


#endif // PUB_HPP
