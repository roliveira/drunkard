#ifndef IO_VALIDATE_H_
#define IO_VALIDATE_H_


#include <vector>
#include <string>


const std::vector<std::string> PARTICLE_ARGS {
    "ic", 
    "bc",
    "n",
    "",
    "",
    "",
    "",
    "",
    ""
};


class ValidateArguments {
private:

    std::vector<std::string> args;

public: 
    
    // Constructor

    ValidateArguments(std::vector<std::string> _args);
    
};

inline ValidateArguments::ValidateArguments(std::vector<std::string> _args) : args(_args) {

}


#endif // IO_VALIDATE_H_
