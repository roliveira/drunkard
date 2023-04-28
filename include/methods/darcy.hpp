
#ifndef METHODS_DARCY_H_
#define METHODS_DARCY_H_


#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <random>
#include <limits>
#include <stdio.h>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include "xtensor/xview.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"

#include "utils/index.hpp"
#include "geometry/geometry.hpp"
#include "methods/flow.hpp"
#include "methods/psi.hpp"


class DarcyMethod : public FlowMethod 
{
public: 

    int ndirections;

    // Viscoity

    double mu;

    // Molecular diffusion

    double dm;

    // Porosity

    xt::xarray<double> phi;
    double phi_mean;
    double phi_stddev;
    double phi_min; 
    double phi_max;

    // Permeability

    xt::xarray<double> perm;
    double perm_mean;
    double perm_stddev;
    double perm_min;
    double perm_max;
    
    // Pressure

    xt::xarray<double> pressure;
    double pressure_mean;
    double pressure_stddev;
    double pressure_min; 
    double pressure_max;

    double pressure_in;
    double pressure_out;

    // Volumetric rate

    double q;

    // Velocity

    xt::xarray<double> velocity;
    xt::xarray<double> velocity_mean;
    xt::xarray<double> velocity_stddev;
    xt::xarray<double> velocity_min;
    xt::xarray<double> velocity_max;

    // Advective time - T1

    xt::xarray<double> t1;
    xt::xarray<double> t1_mean;
    xt::xarray<double> t1_stddev;
    xt::xarray<double> t1_min;
    xt::xarray<double> t1_max;

    // Diffusive time - T2

    xt::xarray<double> t2;
    xt::xarray<double> t2_mean;
    xt::xarray<double> t2_stddev;
    xt::xarray<double> t2_min;
    xt::xarray<double> t2_max;

    // Heterogeneity factor - beta

    xt::xarray<double> beta;
    double beta_mean;
    double beta_stddev;
    double beta_min;
    double beta_max;

    // Pe

    xt::xarray<double> pe;
    xt::xarray<double> pe_mean;
    xt::xarray<double> pe_stddev;
    xt::xarray<double> pe_min;
    xt::xarray<double> pe_max;

    // Direction probability tensor - pij

    xt::xarray<double> direction_prob;
    xt::xarray<double> direction_prob_cumm;

    // Transit Time distribution - psi

    xt::xarray<Psi> psi;
    xt::xarray<PsiADE> psi_ade;

    // Constructor

    DarcyMethod(void);
    DarcyMethod(std::vector<int> shape, std::vector<double> size, std::vector<double> origin);

    // Setters

    void SetBeta(double beta);
    void SetPressure(double _pressure_in, double _pressure_out);

    void SetPorosity(std::vector<double> _phi);
    void SetPorosity(double _phi, int i, int j, int k);

    void SetPermeability(std::vector<double> _perm);
    void SetPermeability(double _perm, int i, int j, int k);

    // Methods

    void CreateNormalPorosity(double  _phi_mean, double  _phi_stddev, int _seed);
    void CreateCorrelatedPorosity(double  _phi_mean, double  _phi_stddev, int _seed);
    void CreateLogNormalPermeability (double _perm_mean, double _perm_stddev, int _seed);

    void CalculateCorrelatedPermeability(double a, double b);
    void CalculateKozenyCarmanPermeability(double sv, double kappa);

    void LabelBoundaryConditions(void);

    void SolvePressure(MPI_Comm comm, int argc, char** argv);
    void CalculateVelocity(void);
    void CalculateProbability(void);
    void CalculateProbabilityADE(void);

    void EnforceVelocity(double vm);
    
    void Summary(void);

    // I/O

    void WritePorosityField(std::string _fname);
    void WritePermeabilityField(std::string _fname);
    void WritePressureField(std::string _fname);
    void WriteVelocityField(std::string _fname);
    void WriteLocationTensor(std::string _fname);
};


// 
// Constructors
//

DarcyMethod::DarcyMethod(void) 
: FlowMethod() 
{}

DarcyMethod::DarcyMethod(std::vector<int> shape, std::vector<double> size, std::vector<double> origin) 
: FlowMethod(shape, size, origin)
{}

//
// Setters
//

void DarcyMethod::SetBeta(double _beta)
{
    beta = xt::empty<double>(this->_shape);
    beta.fill(_beta);
}

void DarcyMethod::SetPressure(double _pressure_in, double _pressure_out) 
{
    pressure_in  = _pressure_in;
    pressure_out = _pressure_out;
}

void DarcyMethod::SetPorosity(std::vector<double> _phi)
{
    phi = xt::adapt(_phi, this->_shape);
}

void DarcyMethod::SetPorosity(double _phi, int i, int j, int k) 
{
    phi(i, j, k) = _phi;
}

void DarcyMethod::SetPermeability(std::vector<double> _perm) 
{
    perm = xt::adapt(_perm, this->_shape);
}

void DarcyMethod::SetPermeability(double _perm, int i, int j, int k) 
{
    perm(i, j, k) = _perm;
}


//
// Methods
//
/*
void DarcyMethod::CreateNormalPorosity(double _phi_mean, double _phi_stddev, int _seed) 
{
    phi_mean   = _phi_mean;
    phi_stddev = _phi_stddev;
    phi        = this->CreateNormalField(phi_mean, phi_stddev, this->_length, _seed);
}

void DarcyMethod::CreateCorrelatedPorosity(double _phi_mean, double _phi_stddev, int _seed) 
{
    phi_mean   = _phi_mean;
    phi_stddev = _phi_stddev;
    phi        = this->CreateCorrelatedField(phi_mean, phi_stddev, this->_length, _seed);
}

void DarcyMethod::CreateLogNormalPermeability(double _perm_mean, double _perm_stddev, int _seed) 
{
    perm_mean   = _perm_mean;
    perm_stddev = _perm_stddev;
    perm        = this->CreateLogNormalField(perm_mean, perm_stddev, this->_length, _seed);
}
*/
void DarcyMethod::CalculateCorrelatedPermeability(double a, double b)
{
    this->perm = a*xt::pow(this->phi, b);
}
/*
void DarcyMethod::CalculateKozenyCarmanPermeability(double sv, double kappa) 
{
    perm    = std::vector<double>(phi.size(), 0.0);
    perm[0] = std::pow(phi[0], 3.0)/(std::pow(sv, 2.0)*kappa*std::pow(1-phi[0], 2.0));

    perm_mean = perm[0];

    for (int i=1; i<perm.size(); ++i) 
    {
        perm[i]   = std::pow(phi[i], 3.0)/(std::pow(sv, 2.0)*kappa*std::pow(1-phi[i], 2.0));
        perm_mean = std::pow(perm_mean, i/(i+1.0)) * std::pow(perm[i], 1.0/(i+1.0));
    }

}
*/
void DarcyMethod::LabelBoundaryConditions(void) 
{
    pressure = xt::zeros<double>(this->_shape);

    for (int i=0; i<this->_length; ++i) 
    {
        std::vector<int> ind = ind2sub(i, this->_shape);

        if (ind[0] == 0)
        {
            this->SetBoundaryCondition(i, PIN);
            pressure(ind[0], ind[1], ind[2]) = pressure_in;
        }
        else if (ind[0] == this->_shape.at(0)-1)
        {
            this->SetBoundaryCondition(i, POUT);
            pressure(ind[0], ind[1], ind[2]) = pressure_out;
        }

    }
}

void DarcyMethod::SolvePressure(MPI_Comm comm, int argc, char** argv) 
{
    // Solve for the pressure
    //     Ax = b,
    // where A is the conductivity matrix, b a vector containing the boundary
    // conditions
    // and x is the pressure.

    Vec vout;  // Temporary vector, to move the results back to the network class

    KSP        ksp;
    PC         pc;
    VecScatter ctx;

    PetscInt istart, istop;

    PetscInt its;
    PetscScalar val;

    PetscErrorCode ierr;

    double minval, maxval, meanval;

    PETSC_COMM_WORLD = comm;

    PetscInitialize(&argc, &argv, NULL, NULL);

    int rank;
    MPI_Comm_rank(comm, &rank);

    int size;
    MPI_Comm_size(comm, &size);

    PetscInt M = this->_length;
    Mat      A;
    
    ierr = MatCreateAIJ(
        comm,
        PETSC_DECIDE, PETSC_DECIDE,
        M, M,
        7, NULL,
        7, NULL,
        &A
    ); // CHKERRQ(ierr);

    ierr = MatSetOption(
        A, 
        MAT_NEW_NONZERO_ALLOCATION_ERR, 
        PETSC_FALSE
    );// CHKERRQ(ierr);

    Vec b;

    ierr = VecCreate(comm, &b);// CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, M);// CHKERRQ(ierr);
    ierr = VecSetType(b, VECMPI); // CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A, &istart, &istop);
    // CHKERRQ(ierr);

    int n = 0;
    std::vector<int> kk;
    std::vector<double> bb;

    for (int i=istart; i<istop; ++i) 
    {
        int m = 0;

        std::vector<int> ii, jj;
        std::vector<double> vv;

        std::vector<int> iind = this->Ind2Sub(i);

        if (this->_flag_bc.at(i) != NONE) 
        {
            ii.push_back(i);
            jj.push_back(i);
            vv.push_back(1);
            ++m;

            kk.push_back(i);
            bb.push_back(this->pressure(iind[0], iind[1], iind[2]));
            ++n;
        }
        else 
        {
            std::vector<int> neighbours = 
                face_neighbours_no_boundary_index(3, i, this->_shape);

            for (int j : neighbours) 
            {
                std::vector<int> jind = this->Ind2Sub(j);

                double temp_perm = std::sqrt(this->perm(iind[0], iind[1], iind[2]))*std::sqrt(this->perm(jind[0], jind[1], jind[2]));
                double g = temp_perm*this->_spacing.at(1)*this->_spacing.at(2)/(this->mu*this->_spacing.at(0));

                ii.push_back(i);
                jj.push_back(j);
                vv.push_back(-g);
                ++m;

                ii.push_back(i);
                jj.push_back(i);
                vv.push_back(g);
                ++m;
            }
        }
        
        ierr = MatSetValues(A, 1, ii.data(), m, jj.data(), vv.data(), ADD_VALUES); // CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);// CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);// CHKERRQ(ierr);
    
    ierr = VecSetValues(b, n, kk.data(), bb.data(), INSERT_VALUES);// CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);// CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);// CHKERRQ(ierr);

    ierr = KSPCreate(comm, &ksp);// CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);// CHKERRQ(ierr);

    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCBJACOBI);

    ierr = KSPSetType(ksp, KSPGMRES);// CHKERRQ(ierr);
    ierr = KSPSetTolerances(
        ksp, 
        std::numeric_limits<double>::epsilon()*std::abs(pressure_in-pressure_out), 
        PETSC_DEFAULT, 
        PETSC_DEFAULT, 
        PETSC_DEFAULT
    );// CHKERRQ(ierr);

    ierr = KSPGetPC(ksp, &pc);//CHKERRQ(ierr);
    ierr = PCSetType(pc, PCSOR);// CHKERRQ(ierr);
    ierr = KSPSolve(ksp, b, b);// CHKERRQ(ierr);

    double norm;
    KSPConvergedReason reason;

    ierr = KSPGetIterationNumber(ksp, &its);// CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &norm);
    ierr = KSPGetConvergedReason(ksp, &reason);

    // std::cout << "iter.: " << its    << std::endl;
    // std::cout << "error: " << norm   << std::endl;
    // std::cout << "reas.: " << reason << std::endl;

    VecScatterCreateToAll(b, &ctx, &vout);
    VecScatterBegin(ctx, b, vout, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, b, vout, INSERT_VALUES, SCATTER_FORWARD);

    for (int i = 0; i < M; ++i) 
    {
        std::vector<int> iind = this->Ind2Sub(i);
        ierr = VecGetValues(vout, 1, &i, &val);// CHKERRQ(ierr);
        pressure(iind[0], iind[1], iind[2]) = val;
    }

    ierr = VecMin (b, NULL,   &minval );// CHKERRQ(ierr);
    ierr = VecMax (b, NULL,   &maxval );// CHKERRQ(ierr);
    ierr = VecNorm(b, NORM_1, &meanval);// CHKERRQ(ierr);

    meanval /= M;

    pressure_min  = minval;
    pressure_max  = maxval;
    pressure_mean = meanval;
    
    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&vout);
    KSPDestroy(&ksp);
    VecScatterDestroy(&ctx);

    PetscFinalize();
}

void DarcyMethod::CalculateVelocity(void) 
{
    std::vector<int> shape = {
        this->_shape.at(0),
        this->_shape.at(1),
        this->_shape.at(2),
        6 
    };

    this->velocity = xt::zeros<double>(shape);
    this->t1 = xt::zeros<double>(shape);
    this->t2 = xt::zeros<double>(shape);
    this->pe = xt::zeros<double>(shape);

    xt::xarray<int> direction = direction_vector();

    for (int k=1; k<this->_shape.at(2)-1; ++k) {
        for (int j=1; j<this->_shape.at(1)-1; ++j) {
            for (int i=1; i<this->_shape.at(0)-1; ++i) {
                xt::xarray<int> iind {i, j, k};

                for (int idir=0; idir<6; ++idir) {
                    xt::xarray<int> jind = xt::view(direction, idir) + iind;

                    double l = this->CalculateDistance(iind, jind);

                    this->t2(i, j, k, idir) = l*l/this->dm;

                    double temp_phi = 0.5*(this->phi(i, j, k)+this->phi(jind[0], jind[1], jind[2]));
                    double temp_perm = std::sqrt(this->perm(i, j, k)*this->perm(jind[0], jind[1], jind[2]));

                    if (temp_phi > 1E-12) {
                        this->velocity(i, j, k, idir) = temp_perm*(this->pressure(jind[0], jind[1], jind[2])-this->pressure(i, j, k))/this->mu/l;
                    }
                    else {
                        this->velocity(i, j, k, idir) = 0;
                    }

                    if (std::abs(velocity(i, j, k, idir)) < 1E-16) {
                        this->t1(i, j, k, idir) = std::numeric_limits<double>::infinity();
                        this->pe(i, j, k, idir) = 0.0;
                    }
                    else {
                        this->t1(i, j, k, idir) = l/std::abs(this->velocity(i, j, k, idir));
                        this->pe(i, j, k, idir) = 2 * this->t2(i, j, k, idir) / this->t1(i, j, k, idir);
                    }
                }
            }   
        }   
    }

    // ghost layer
    
    for (int k=0; k<this->_shape.at(2)-1; ++k) {
        for (int j=0; j<this->_shape.at(1)-1; ++j) {
            for (int i: {0, this->_shape.at(0)-1}) {
                std::vector<int> iind = {i, j, k};
                int ir = this->GetReflectedGridIndex(iind);
                std::vector<int> rind = this->Ind2Sub(ir);

                for (int idir=0; idir<6; ++idir)
                {
                    this->velocity(iind[0], iind[1], iind[2], idir) = this->velocity(rind[0], rind[1], rind[2], idir);
                    this->t1(iind[0], iind[1], iind[2], idir) = this->t1(rind[0], rind[1], rind[2], idir);
                    this->t2(iind[0], iind[1], iind[2], idir) = this->t2(rind[0], rind[1], rind[2], idir);
                    this->pe(iind[0], iind[1], iind[2], idir) = this->pe(rind[0], rind[1], rind[2], idir);
                }
            }
        }
    }

    for (int k=0; k<this->_shape.at(2)-1; ++k) {
        for (int i=0; i<this->_shape.at(0)-1; ++i) {
            for (int j: {0, this->_shape.at(1)-1}) {
                std::vector<int> iind = {i, j, k};
                int ir = this->GetReflectedGridIndex(iind);
                std::vector<int> rind = this->Ind2Sub(ir);

                for (int idir=0; idir<6; ++idir)
                {
                    this->velocity(iind[0], iind[1], iind[2], idir) = this->velocity(rind[0], rind[1], rind[2], idir);
                    this->t1(iind[0], iind[1], iind[2], idir) = this->t1(rind[0], rind[1], rind[2], idir);
                    this->t2(iind[0], iind[1], iind[2], idir) = this->t2(rind[0], rind[1], rind[2], idir);
                    this->pe(iind[0], iind[1], iind[2], idir) = this->pe(rind[0], rind[1], rind[2], idir);
                }
            }
        }
    }

    for (int j=0; j<this->_shape.at(1)-1; ++j) {
        for (int i=0; i<this->_shape.at(0)-1; ++i) {
            for (int k: {0, this->_shape.at(2)-1}) {
                std::vector<int> iind = {i, j, k};
                int ir = this->GetReflectedGridIndex(iind);
                std::vector<int> rind = this->Ind2Sub(ir);

                for (int idir=0; idir<6; ++idir)
                {
                    this->velocity(iind[0], iind[1], iind[2], idir) = this->velocity(rind[0], rind[1], rind[2], idir);
                    this->t1(iind[0], iind[1], iind[2], idir) = this->t1(rind[0], rind[1], rind[2], idir);
                    this->t2(iind[0], iind[1], iind[2], idir) = this->t2(rind[0], rind[1], rind[2], idir);
                    this->pe(iind[0], iind[1], iind[2], idir) = this->pe(rind[0], rind[1], rind[2], idir);
                }
            }
        }
    }

    // for (int i=0; i<this->_length; ++i) 
    // {
    //     if (!this->IsBoundary(i)) continue;

    //     int ir = this->GetReflectedGridIndex(i);

    //     std::vector<int> iind = this->Ind2Sub(i);
    //     std::vector<int> rind = this->Ind2Sub(ir);

    //     for (int idir=0; idir<6; ++idir)
    //     {
    //         this->velocity(iind[0], iind[1], iind[2], idir) = this->velocity(rind[0], rind[1], rind[2], idir);
    //         this->t1(iind[0], iind[1], iind[2], idir) = this->t1(rind[0], rind[1], rind[2], idir);
    //         this->t2(iind[0], iind[1], iind[2], idir) = this->t2(rind[0], rind[1], rind[2], idir);
    //         this->pe(iind[0], iind[1], iind[2], idir) = this->pe(rind[0], rind[1], rind[2], idir);
    //     }
    // }
}

void DarcyMethod::CalculateProbability(void)
{
    float eps         = std::numeric_limits<float>::epsilon();
    int dimension     = 3;
    this->ndirections = 6;
    int nsize         = velocity.size() / this->ndirections;
    
    std::vector<int> shape = {
        this->_shape.at(0),
        this->_shape.at(1),
        this->_shape.at(2),
        6 
    };

    this->psi = xt::empty<Psi>(shape);

    xt::xarray<int> direction = direction_vector();

    this->direction_prob      = xt::zeros<double>(shape);
    this->direction_prob_cumm = xt::zeros<double>(shape);

    for (int k=1; k<this->_shape.at(2)-1; ++k)
    {
        for (int j=1; j<this->_shape.at(1)-1; ++j)
        {
            for (int i=1; i<this->_shape.at(0)-1; ++i)
            {
                xt::xarray<int> iind {i, j, k};

                for (int idir=0; idir<6; ++idir)
                {
                    xt::xarray<int> jind = xt::view(direction, idir) + iind;

                    this->psi(i, j, k, idir) = 
                        Psi(this->t1(i, j, k, idir), this->t2(i, j, k, idir), this->beta(i, j, k));

                    // pe ~~ 0 : pe / (1 - e^-pe) -> 1,
                    //           pe / (e^+pe - 1) -> 1.
                    if (std::abs(this->pe(i, j, k, idir)) < eps * 1.0)
                    {
                        this->direction_prob(i, j, k, idir) = 1.0;
                    }
                    // pe >> 1 : pe / (1 - e^-pe) -> pe,
                    //           pe / (e^+pe - 1) -> 0.
                    else if (this->pe(i, j, k, idir) > 50) 
                    {
                        if (this->velocity(i, j, k, idir) > 0) 
                            this->direction_prob(i, j, k, idir) = this->pe(i, j, k, idir);
                        else                                 
                            this->direction_prob(i, j, k, idir) = 0.0;
                    }
                    else
                    {
                        if (this->velocity(i, j, k, idir) > 0)
                        {
                            this->direction_prob(i, j, k, idir) = 
                                this->pe(i, j, k, idir) / (1 - std::exp(-this->pe(i, j, k, idir)));
                        }
                        else
                        {
                            this->direction_prob(i, j, k, idir) = 
                                this->pe(i, j, k, idir) / (std::exp(+this->pe(i, j, k, idir)) - 1);
                        }
                    }
                }
            }
        }
    }

    for (int k=1; k<this->_shape.at(2)-1; ++k)
    {
        for (int j=1; j<this->_shape.at(1)-1; ++j)
        {
            for (int i=1; i<this->_shape.at(0)-1; ++i)
            {
                this->direction_prob_cumm(i, j, k, 0) = this->direction_prob(i, j, k, 0);

                for (int idir=1; idir<6; ++idir)
                {
                    this->direction_prob_cumm(i, j, k, idir) = this->direction_prob_cumm(i, j, k, idir-1) + this->direction_prob(i, j, k, idir);
                }

                for (int idir=0; idir<6; ++idir)
                {
                    this->direction_prob(i, j, k, idir)      /= this->direction_prob_cumm(i, j, k, 5);
                    this->direction_prob_cumm(i, j, k, idir) /= this->direction_prob_cumm(i, j, k, 5);
                }
            }
        }
    }

    for (int i = 0; i < nsize; ++i) 
    {
        if (!this->IsBoundary(i)) continue;

        int ir = this->GetReflectedGridIndex(i);

        std::vector<int> iind = this->Ind2Sub(i);
        std::vector<int> rind = this->Ind2Sub(ir);

        for (int idir=0; idir<6; ++idir) 
        {
            this->psi(iind[0], iind[1], iind[2], idir) = 
                this->psi(rind[0], rind[1], rind[2], idir);

            this->direction_prob(iind[0], iind[1], iind[2], idir) = 
                this->direction_prob(rind[0], rind[1], rind[2], idir);

            this->direction_prob_cumm(iind[0], iind[1], iind[2], idir) = 
                this->direction_prob_cumm(rind[0], rind[1], rind[2], idir);
        }
    }
    
}

void DarcyMethod::CalculateProbabilityADE(void)
{
    float eps     = std::numeric_limits<float>::epsilon();
    int dimension = 3;
    ndirections   = 6;
    int nsize     = velocity.size() / ndirections;
    
    std::vector<int> shape = {
        this->_shape.at(0),
        this->_shape.at(1),
        this->_shape.at(2),
        6 
    };

    this->psi_ade = xt::empty<PsiADE>(this->_shape);

    xt::xarray<int> direction = direction_vector();

    direction_prob = xt::zeros<double>(shape);
    direction_prob_cumm = xt::zeros<double>(shape);

    for (int k=1; k<this->_shape.at(2)-1; ++k)
    {
        for (int j=1; j<this->_shape.at(1)-1; ++j)
        {
            for (int i=1; i<this->_shape.at(0)-1; ++i)
            {
                xt::xarray<int> iind {i, j, k};

                std::vector<double> temp_v(6, 0.0);
                std::vector<double> temp_pe(6, 0.0);
                std::vector<double> temp_l(6, 0.0);

                for (int idir=0; idir<6; ++idir)
                {
                    xt::xarray<int> jind = xt::view(direction, idir) + iind;

                    temp_v[idir] = velocity(i, j, k, idir);
                    temp_pe[idir] = pe(i, j, k, idir);

                    temp_l[idir] = this->GetDistance(
                        std::vector<int>(iind.begin(), iind.end()),
                        std::vector<int>(jind.begin(), jind.end())
                    );

                    // pe ~~ 0 : pe / (1 - e^-pe) -> 1,
                    //           pe / (e^+pe - 1) -> 1.
                    if (std::abs(pe(i, j, k, idir)) < eps * 1.0)
                    {
                        direction_prob(i, j, k, idir) = 1.0;
                    }
                    // pe >> 1 : pe / (1 - e^-pe) -> pe,
                    //           pe / (e^+pe - 1) -> 0.
                    else if (pe(i, j, k, idir) > 50) 
                    {
                        if (velocity(i, j, k, idir) > 0) 
                            direction_prob(i, j, k, idir) = pe(i, j, k, idir);
                        else                                 
                            direction_prob(i, j, k, idir) = 0.0;
                    }
                    else 
                    {
                        if (velocity(i, j, k, idir) > 0) 
                            direction_prob(i, j, k, idir) = pe(i, j, k, idir) / (1 - std::exp(-pe(i, j, k, idir)));
                        else                                 
                            direction_prob(i, j, k, idir) = pe(i, j, k, idir) / (std::exp(+pe(i, j, k, idir)) - 1);
                    }
                }

                psi_ade(i, j, k) = PsiADE(temp_v, this->dm, temp_pe, temp_l);
            }
        }
    }

    for (int k=1; k<this->_shape.at(2)-1; ++k)
    {
        for (int j=1; j<this->_shape.at(1)-1; ++j)
        {
            for (int i=1; i<this->_shape.at(0)-1; ++i)
            {
                direction_prob_cumm(i, j, k, 0) = direction_prob(i, j, k, 0);

                for (int idir=1; idir<6; ++idir)
                {
                    direction_prob_cumm(i, j, k, idir) = direction_prob_cumm(i, j, k, idir-1) + direction_prob(i, j, k, idir);
                }

                for (int idir=0; idir<6; ++idir)
                {
                    direction_prob_cumm(i, j, k, idir) /= direction_prob_cumm(i, j, k, 5);
                }
            }
        }
    }

    for (int i = 0; i < nsize; ++i) 
    {
        if (!this->IsBoundary(i)) continue;

        int ir = this->GetReflectedGridIndex(i);

        std::vector<int> iind = this->Ind2Sub(i);
        std::vector<int> rind = this->Ind2Sub(ir);

        psi_ade(iind[0], iind[1], iind[2]) = psi_ade(rind[0], rind[1], rind[2]);

        for (int idir=0; idir<6; ++idir) 
        {
            direction_prob(iind[0], iind[1], iind[2], idir) = direction_prob(rind[0], rind[1], rind[2], idir);
            direction_prob_cumm(iind[0], iind[1], iind[2], idir) = direction_prob_cumm(rind[0], rind[1], rind[2], idir);
        }
    }
    
}

void DarcyMethod::Summary(void) {

    // Basic stats

    this->phi_mean = 0.0;
    this->phi_min  = std::numeric_limits<double>::max();
    this->phi_max  = std::numeric_limits<double>::lowest();

    this->perm_mean = 0.0;
    this->perm_min  = std::numeric_limits<double>::max();
    this->perm_max  = std::numeric_limits<double>::lowest();

    this->pressure_mean = 0.0;
    this->pressure_min  = std::numeric_limits<double>::max();
    this->pressure_max  = 0.0;

    this->velocity_min  = xt::zeros<double>({3});
    this->velocity_max  = xt::zeros<double>({3});
    this->velocity_min += std::numeric_limits<double>::max();
    this->velocity_max += std::numeric_limits<double>::lowest();

    this->t1_min  = xt::zeros<double>({3});
    this->t1_max  = xt::zeros<double>({3});
    this->t1_min += std::numeric_limits<double>::max();
    this->t1_max += std::numeric_limits<double>::lowest();

    this->t2_min  = xt::zeros<double>({3});
    this->t2_max  = xt::zeros<double>({3});
    this->t2_min += std::numeric_limits<double>::max();
    this->t2_max += std::numeric_limits<double>::lowest();

    this->pe_min  = xt::zeros<double>({3});
    this->pe_max  = xt::zeros<double>({3});
    this->pe_min += std::numeric_limits<double>::max();
    this->pe_max += std::numeric_limits<double>::lowest();

    for (int i=0; i<this->_length; ++i) {
        std::vector<int> iind = this->Ind2Sub(i);

        double temp_phi   = this->phi(iind[0], iind[1], iind[2]);
        double temp_perm  = this->perm(iind[0], iind[1], iind[2]);
        double temp_press = this->pressure(iind[0], iind[1], iind[2]);

        this->phi_mean      = this->phi_mean+(temp_phi-this->phi_mean)/(i+1);
        this->perm_mean     = this->perm_mean+(temp_perm-this->perm_mean)/(i+1);
        this->pressure_mean = this->pressure_mean+(temp_press-this->pressure_mean)/(i+1);

        if (temp_phi   < this->phi_min)      this->phi_min      = temp_phi;
        if (temp_phi   > this->phi_max)      this->phi_max      = temp_phi;
        if (temp_perm  < this->perm_min)     this->perm_min     = temp_perm;
        if (temp_perm  > this->perm_max)     this->perm_max     = temp_perm;
        if (temp_press < this->pressure_min) this->pressure_min = temp_press;
        if (temp_press > this->pressure_max) this->pressure_max = temp_press;

        for (int d=0; d<3; ++d) {
            double temp_v  = this->velocity(iind[0], iind[1], iind[2], 2*d+1);
            double temp_t1 = this->t1(iind[0], iind[1], iind[2], 2*d+1);
            double temp_t2 = this->t2(iind[0], iind[1], iind[2], 2*d+1);
            double temp_pe = this->pe(iind[0], iind[1], iind[2], 2*d+1);

            if (temp_v  < this->velocity_min(d)) this->velocity_min(d) = temp_v;
            if (temp_v  > this->velocity_max(d)) this->velocity_max(d) = temp_v;
            if (temp_t1 < this->t1_min(d))       this->t1_min(d)       = temp_t1;
            if (temp_t1 > this->t1_max(d))       this->t1_max(d)       = temp_t1;
            if (temp_t2 < this->t2_min(d))       this->t2_min(d)       = temp_t2;
            if (temp_t2 > this->t2_max(d))       this->t2_max(d)       = temp_t2;
            if (temp_pe < this->pe_min(d))       this->pe_min(d)       = temp_pe;
            if (temp_pe > this->pe_max(d))       this->pe_max(d)       = temp_pe;
        }
    }

    velocity_mean    = xt::zeros<double>({3});
    velocity_mean(0) = (perm_mean/mu)*(pressure_out-pressure_in)/this->_size.at(0);
    t1_mean          = xt::zeros<double>({3});
    t1_mean(0)       = this->_size.at(0)/std::abs(velocity_mean[0]);
    t2_mean          = xt::zeros<double>({3});
    t2_mean(0)       = this->_size.at(0)*this->_size.at(0)/std::abs(dm);
    pe_mean          = xt::zeros<double>({3});
    pe_mean(0)       = 2 * t2_mean(0) / t1_mean(0);

    this->q = 0;

    for (int k=0; k<this->_shape[2]; ++k)
    {
        for (int j=0; j<this->_shape[1]; ++j)
        {
            this->q += this->_spacing[1]*this->_spacing[2]*this->velocity(0, j, k, 1);
        }
    }

    // this->q = velocity_mean(0)*this->_size.at(1)*this->_size.at(2);
}

void DarcyMethod::EnforceVelocity(double vm)
{
    std::vector<int> shape = {
        this->_shape.at(0),
        this->_shape.at(1),
        this->_shape.at(2),
        6
    };

    this->velocity = xt::zeros<double>(shape);
    this->t1 = xt::zeros<double>(shape);
    this->t2 = xt::zeros<double>(shape);
    this->pe = xt::zeros<double>(shape);

    // (-1, 0, 0), (+1, 0, 0), (0, -1, 0), ...
    std::vector<double> v(6, 0);
    v[0] = -vm;
    v[1] = +vm;

    xt::xarray<int> direction = direction_vector();

    for (int k=0; k<shape[2]; ++k)
    {
        for (int j=0; j<shape[1]; ++j)
        {
            for (int i=0; i<shape[0]; ++i)
            {
                xt::xarray<int> iind {i, j, k};

                for (int d=0; d<shape[3]; ++d)
                {
                    xt::xarray<int> jind = xt::view(direction, d) + iind;

                    double l = this->GetDistance(
                        std::vector<int>(iind.begin(), iind.end()), 
                        std::vector<int>(jind.begin(), jind.end())
                    );

                    this->velocity(i, j, k, d) = v[d];

                    if (std::abs(this->velocity(i, j, k, d)) < 1E-12)
                    {
                        this->velocity(i, j, k, d) = 0.0;
                        this->t1(i, j, k, d) = std::numeric_limits<double>::infinity();
                        this->t2(i, j, k, d) = l*l/this->dm;
                        this->pe(i, j, k, d) = 0.0;
                    }
                    else
                    {
                        this->t1(i, j, k, d) = l/std::abs(this->velocity(i, j, k, d));
                        this->t2(i, j, k, d) = l*l/this->dm;
                        this->pe(i, j, k, d) = 2.0*this->t2(i, j, k, d)/this->t1(i, j, k, d);
                    }

                    std::cout << "v  = " << this->velocity(i, j, k, d) << std::endl;
                    std::cout << "t1 = " << this->t1(i, j, k, d) << std::endl;
                    std::cout << "t2 = " << this->t2(i, j, k, d) << std::endl;
                    std::cout << "pe = " << this->pe(i, j, k, d) << std::endl;
                }
            }
        }
    }
}

// I/O

void DarcyMethod::WritePorosityField(std::string _fname) 
{
    ToVTI(
        _fname.c_str(), 
        this->phi, 
        this->_shape, 
        this->_spacing, 
        this->_origin
    );
}

void DarcyMethod::WritePermeabilityField (std::string _fname) 
{
    ToVTI(
        _fname.c_str(), 
        this->perm, 
        this->_shape, 
        this->_spacing, 
        this->_origin
    );
}

void DarcyMethod::WritePressureField(std::string _fname) 
{
    ToVTI(
        _fname.c_str(), 
        this->pressure, 
        this->_shape, 
        this->_spacing, 
        this->_origin
    );
}

void DarcyMethod::WriteLocationTensor(std::string _fname) 
{
    std::vector<double> temp = std::vector<double>(direction_prob.begin(), direction_prob.end());

    // Tensor2VTI(
    //     (_fname + "_loc_prob.vti").c_str(),
    //     this->direction_prob,
    //     this->_shape,
    //     this->_spacing,
    //     this->_origin
    // );

    temp = std::vector<double>(direction_prob_cumm.begin(), direction_prob_cumm.end());

    Tensor2VTI(
        // (_fname + "_loc_prob_cumm.vti").c_str(),
        _fname.c_str(),
        this->direction_prob_cumm,
        this->_shape,
        this->_spacing,
        this->_origin
    );

    // Tensor2VTI(
    //     (_fname + "_loc_prob.vti").c_str(),
    //     this->direction_prob,
    //     this->_shape,
    //     this->_spacing,
    //     this->_origin
    // );

    // Tensor2VTI(
    //     (_fname + "_loc_prob_cumm.vti").c_str(),
    //     this->direction_prob_cumm,
    //     this->_shape,
    //     this->_spacing,
    //     this->_origin
    // );
}


std::ostream &operator<<(std::ostream &os, const std::shared_ptr<DarcyMethod> &tm) 
{
    os << static_cast<const std::shared_ptr<Geometry>>(tm);

    os << message("Viscosity [Pa.s]", tm->mu);

    os << message("Poro. mean [-]", tm->phi_mean);
    os << message("Poro. min. [-]", tm->phi_min);
    os << message("Poro. max. [-]", tm->phi_max);

    os << message("Perm. mean [m2]", tm->perm_mean);
    os << message("Perm. min. [m2]", tm->perm_min);
    os << message("Perm. max. [m2]", tm->perm_max);

    os << message("Pressure mean [Pa]", tm->pressure_mean);
    os << message("Pressure min. [Pa]", tm->pressure_min);
    os << message("Pressure max. [Pa]", tm->pressure_max);

    os << message("Mean flux [m3/s]", tm->q);

    os << message("Velocity mean [m/s]", std::vector<double>(tm->velocity_mean.begin(), tm->velocity_mean.end()));
    os << message("Velocity min. [m/s]", std::vector<double>(tm->velocity_min.begin(), tm->velocity_min.end()));
    os << message("Velocity max. [m/s]", std::vector<double>(tm->velocity_max.begin(), tm->velocity_max.end()));
    
    os << message("Adv. transit time mean [s]", std::vector<double>(tm->t1_mean.begin(), tm->t1_mean.end()));
    os << message("Adv. transit time min. [s]", std::vector<double>(tm->t1_min.begin(), tm->t1_min.end()));
    os << message("Adv. transit time max. [s]", std::vector<double>(tm->t1_max.begin(), tm->t1_max.end()));
    
    os << message("Cut-off diff. time mean [s]", std::vector<double>(tm->t2_mean.begin(), tm->t2_mean.end()));
    os << message("Cut-off diff. time min. [s]", std::vector<double>(tm->t2_min.begin(), tm->t2_min.end()));
    os << message("Cut-off diff. time max. [s]", std::vector<double>(tm->t2_max.begin(), tm->t2_max.end()));

    return os;
}



#endif  // METHODS_DARCY_H_
