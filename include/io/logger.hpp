#ifndef IO_LOGGER_H_
#define IO_LOGGER_H_


#include <cstdio>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <mpi.h>

#include "io/echo.hpp"


class Logger 
{

private:

    MPI_Comm _comm;

    const char    *fname;
    std::ofstream fout; 

    std::clock_t begin;

public:

    // Constructor

    Logger(const char *fname, MPI_Comm comm);

    // Destructor

    ~Logger(void);

    // Stream operator

    template <class T> 
    Logger &operator<<(const T &msg);
    Logger &operator<<(std::ostream& (*os)(std::ostream&));

    // Methods

    double GetElapsedTime(void);
};

inline Logger::Logger(const char *fname, MPI_Comm comm)
: begin(std::clock()), _comm(comm)
{
    int rank;
    MPI_Comm_rank(_comm, &rank);

    if (rank == 0)
    {
        fout = std::ofstream(fname);

        std::cout << banner() << std::endl;
        fout      << banner() << std::endl;
    }
}

inline Logger::~Logger(void) 
{
    int rank;
    MPI_Comm_rank(_comm, &rank);
    
    if (rank == 0)
    {
        double elapsed_secs = GetElapsedTime();
        std::string out = section_close("Simulation completed", elapsed_secs);

        std::cout << out << std::endl << std::endl;
        fout      << out << std::endl << std::endl;

        fout.close();
    }
}

template <class T> 
inline Logger &Logger::operator<<(const T &msg) 
{
    int rank;
    MPI_Comm_rank(_comm, &rank);
    
    if (rank == 0)
    {
        std::cout << msg;
        fout      << msg;
    }

    return *this;
}

inline Logger &Logger::operator<<(std::ostream& (*os)(std::ostream&)) 
{
    int rank;
    MPI_Comm_rank(_comm, &rank);
    
    if (rank == 0)
    {            
        std::cout << os;
        fout      << os;
    }

    return *this;
}

double Logger::GetElapsedTime(void) 
{
    std::clock_t end = std::clock();
    return double(end - begin) / CLOCKS_PER_SEC;
}


#endif // IO_LOGGER_H_
