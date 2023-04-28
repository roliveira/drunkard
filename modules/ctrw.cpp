
#include <iostream>
#include <memory>

#include <mpi.h>

// #include "nlohmann/json.hpp"

#include "xtensor/xnpy.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xreducer.hpp"

#include "io/echo.hpp"
#include "utils/pub.hpp"
#include "geometry/geometry.hpp"
#include "methods/darcy.hpp"
#include "methods/psi.hpp"


// using json = nlohmann::json;
#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Common.hpp>
using namespace Reaktoro;


int main(int argc, char **argv) 
{

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (argc < 2 && world_rank==0) {    
        std::cerr << "Usage: " << argv[0] << " FILENAME.case" << std::endl;
        return EXIT_FAILURE;
    }

    if (std::string(argv[1]).find(FILE_EXTENSION) == std::string::npos) {
        std::cerr << "Unknown extension " << argv[1] << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream ifs(argv[1]);
    json parameters;
    ifs >> parameters;

    // File names

    std::string fname;

    if (parameters.count("fname")) {
        fname = parameters["fname"];
    }
    else {
        std::size_t pos_ext = std::string(argv[1]).rfind(FILE_EXTENSION);
        fname = std::string(argv[1]).substr(pos_ext);
    }

    std::string lname = std::string(argv[1]) + ".log";

    std::shared_ptr<Pub> pub = std::make_shared<Pub>(comm);
    Logger log((lname).c_str(), comm);

    pub->_fprefix = fname;
    pub->_log = &log;

    log << section("Beginning of the input file");
    log << parameters.dump(4) << std::endl << std::endl;
    log << section("End of the input file");

    log << section("Input parameters");

    log << message("input file name", argv[1]);
    log << message("log file name", lname);
    log << message("output file prefix", fname);

    // Output options

    bool write_vtk = parameters.value("vtk", true);
    bool write_npy = parameters.value("npy", true);

    log << message("write vtk files", write_vtk);
    log << message("write npy files", write_npy);

    // Geometry

    std::vector<int> n = parameters["geometry"]["shape"].get<std::vector<int>>();
    std::vector<double> l = parameters["geometry"]["size"].get<std::vector<double>>();
    std::vector<double> o = parameters["geometry"].value("origin", std::vector<double>(3, 0.0));

    log << message("geometry shape", n);
    log << message("geometry size", l);
    log << message("geometry origin", o);

    bool solve_flow = false;
    std::string fporosity;
    std::string fpermeability;
    std::string fpressure;
    std::vector<double> pressure_bc;

    if (parameters["flow"]["pressure"].count("file")) {
        solve_flow = false;
        fpressure = parameters["flow"]["pressure"]["file"].get<std::string>();
        
        log << message("pressure file", fpressure);
    }
    else {
        solve_flow = true;
        pressure_bc = parameters["flow"]["pressure"]["boundary"].get<std::vector<double>>();

        log << message("pressure bc", pressure_bc);
    }

    double mu = parameters["flow"]["viscosity"].get<double>();
    fporosity = parameters["flow"]["porosity"]["file"].get<std::string>();
    fpermeability = parameters["flow"]["permeability"]["file"].get<std::string>();

    log << message("porosity file", fporosity);
    log << message("permeability file", fpermeability);
    
    log << message("dynamic viscosity", mu);
    
    // Transport solver

    double beta = parameters["transport"]["beta"].get<double>();
    double dm = parameters["transport"]["dm"].get<double>();
    double dt = parameters["transport"]["dt"].get<double>();
    double tmax = parameters["transport"]["tmax"].get<double>();

    int step = parameters["transport"].value("step", 0);
    double t0 = parameters["transport"].value("t0", 0.0);

    log << message("ctrw beta", beta);
    log << message("dm", dm);
    log << message("dt", dt);
    log << message("tmax", tmax);
    log << message("step", step);
    log << message("t0", t0);

    int num_parts = parameters["transport"]["initial"]["n"].get<int>();

    std::string transport_ic_type = parameters["transport"]["initial"]["type"].get<std::string>();
    INITIAL_CONDITION transport_ic = map_initial_condition[transport_ic_type];
    std::vector<int> transport_index = parameters["transport"]["initial"]["index"].get<std::vector<int>>();

    std::string transport_bc_type = parameters["transport"]["boundary"].get<std::string>();
    BOUNDARY_CONDITION transport_bc = map_boundary_condition[transport_bc_type];

    log << message("initial particle number", num_parts);
    log << message("initial particle source", transport_index);
    log << message("particle initial condition", transport_ic_type);
    log << message("particle boundary condition", transport_bc_type);
    
    
    /*
     * Flow solver
     */


    log << section("Flow solver");

    std::shared_ptr<DarcyMethod> darcy = std::make_shared<DarcyMethod>(n, l, o);

    darcy->phi  = xt::load_npy<double>(fporosity.c_str());
    darcy->perm = xt::load_npy<double>(fpermeability.c_str());

    darcy->mu = mu;
    darcy->dm = dm;
    darcy->beta = xt::zeros<double>(darcy->_shape) + beta;

    if (solve_flow) {
        darcy->SetPressure(pressure_bc.at(0), pressure_bc.at(1));
        darcy->LabelBoundaryConditions();
        darcy->SolvePressure(MPI_COMM_WORLD, argc, argv);
    }
    else {
        darcy->pressure = xt::load_npy<double>(fpressure);
    }

    darcy->Summary();
    darcy->CalculateVelocity();
    darcy->CalculateProbability();

    if (write_vtk) {
        darcy->WritePorosityField(    (fname + "_poro0.vti").c_str());
        darcy->WritePermeabilityField((fname + "_perm0.vti").c_str());
        darcy->WritePressureField(    (fname + "_pres0.vti").c_str());
        // darcy->WriteLocationTensor(   (fname + "_dirp0.vti").c_str());
        // darcy->WriteVelocityField(    (fname + "_velo0.vti").c_str());
    }

    if (write_npy) {
        xt::dump_npy((fname+"_poro0.npy").c_str(), darcy->phi                );
        xt::dump_npy((fname+"_perm0.npy").c_str(), darcy->perm               );
        xt::dump_npy((fname+"_pres0.npy").c_str(), darcy->pressure           );
        xt::dump_npy((fname+"_dirp0.npy").c_str(), darcy->direction_prob     );
        xt::dump_npy((fname+"_dirc0.npy").c_str(), darcy->direction_prob_cumm);
        xt::dump_npy((fname+"_velo0.npy").c_str(), darcy->velocity           );
        xt::dump_npy((fname+"_t10.npy"  ).c_str(), darcy->t1                 );
        xt::dump_npy((fname+"_t20.npy"  ).c_str(), darcy->t2                 );
        xt::dump_npy((fname+"_pe0.npy"  ).c_str(), darcy->pe                 );
    }

    log << darcy << std::endl;
    pub->_darcy = std::move(darcy);

    
    /*
     * Transport solver
     */


    log << section("Transport solver");

    if (transport_ic == INITIAL_CONDITION::POINT) {
        pub->CreateParticlePointInjection(num_parts, transport_index[0], 0);
    }
    else if (transport_ic == INITIAL_CONDITION::PLANE) {
       pub->CreateParticlePlaneInjection(num_parts, transport_index[0], 0);
    }
    else if (transport_ic == INITIAL_CONDITION::VOLUME) {
       pub->CreateParticleDistribution(num_parts, 0.0);
    }
    else /* (transport_ic == INITIAL_CONDITION::FLOW) */ {
        pub->CreateParticlePlaneInjectionPhi0(num_parts, transport_index[0], 0.0);
    }

    pub->SetTime(step, t0, dt, tmax);
    pub->SetInitialConditions();
    pub->WriteParticles();
    pub->CreateParticleSummaryFile();


    while (pub->_t < pub->_tmax) {
        for (std::list<Particle*> &plist : pub->_cells._parts) {
            std::list<Particle*>::iterator p = plist.begin();

            while (p != plist.end() && (*p)->_valid) {
                while ((*p)->_dtt < pub->_dt && (*p)->_valid) {
                    xt::xarray<int> source_ind = (*p)->GetSourceI();
                    xt::xarray<int> target_ind = (*p)->GetTargetI();
                    xt::xarray<int> target_dir = xt::view(direction_vector(), (*p)->_dir, xt::all());

                    // 
                    // apply BC
                    // 

                    if (pub->_darcy->IsBoundary(target_ind)) {
                        if (transport_bc == BOUNDARY_CONDITION::PERIODIC) {
                            target_ind = pub->_darcy->GetReflectedGridSubIndex(target_ind);
                        }
                        else /* (transport_bc == BOUNDARY_CONDITION::BT) */ {
                            if (pub->_darcy->IsXp(target_ind)) {
                                ++(pub->_bt);
                                (*p)->_valid = false;
                            }
                            else {
                                target_ind = pub->_darcy->GetBouncedGridSubindex(target_ind);
                            }
                        }
                    }
                    
                    if ((*p)->_valid) {
                        
                        if (transport_bc == BOUNDARY_CONDITION::PERIODIC) {
                            (*p)->dx  = pub->_darcy->GetLocation(target_ind) - pub->_darcy->GetLocation(source_ind);
                            (*p)->xg += (*p)->dx;
                        }
                        else {
                            (*p)->dx = 0.0;
                            (*p)->xg = pub->_darcy->GetLocation(target_ind);
                        }

                        (*p)->SetX((*p)->xg(0), (*p)->xg(1), (*p)->xg(2));
                        (*p)->SetI(target_ind(0), target_ind(1), target_ind(2));

                        (*p)->_dir = pub->SampleDirection(target_ind(0), target_ind(1), target_ind(2));
                        (*p)->_tt = pub->SampleTransitTime(target_ind(0), target_ind(1), target_ind(2), (*p)->_dir);

                        (*p)->_dtt = (*p)->_dtt + (*p)->_tt;
                        (*p)->v = pub->_darcy->CalculateDistance((*p)->GetSourceI(), (*p)->GetTargetI())/(*p)->_tt;

                        (*p)->_move_list = true;
                    }
                }

                if ((*p)->_valid) {
                    xt::xarray<int> target_dir = xt::view(direction_vector(), (*p)->_dir, xt::all());
                    double target_time = pub->_dt - ((*p)->_dtt - (*p)->_tt);

                    (*p)->_t  += pub->_dt;
                    (*p)->_dtt = (*p)->_dtt - pub->_dt;

                    // mid-flight broken!

                    // (*p)->dx  = (*p)->v*target_dir*target_time;
                    // (*p)->xg += (*p)->dx;
                    (*p)->SetX((*p)->xg(0), (*p)->xg(1), (*p)->xg(2));

                    (*p)->_move_list = false;
                }

                if (!(*p)->_valid) {
                    p = plist.erase(p);
                }
                else if ((*p)->_move_list) {
                    (*p)->_move_list = false;

                    std::list<Particle*> &plist_target = pub->_cells._parts((*p)->_i,  (*p)->_j,  (*p)->_k );
                    (*p)->SetI0((*p)->_i, (*p)->_j, (*p)->_k);
                    
                    plist_target.splice(plist_target.begin(), plist, p++);
                }
                else {
                    ++p;
                }
            }
        }

        pub->_step += 1;
        pub->_t += dt;

        pub->CalculateParticleSummary();
        pub->WriteParticleSummary();
        pub->WriteParticles();

        *(pub->_log) << pub->OutputMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
