#ifndef IO_PARSER_H_
#define IO_PARSER_H_


#include <string>
#include <vector>
#include <map>
#include <exception>
#include <stdint.h>
#include <iostream> 
#include <utility>

#include "cpptoml.h"


class GeneralParameters {
private:
    std::shared_ptr<cpptoml::table> table;
    
public:
    std::string fname;

    // Constructor

    GeneralParameters(const char *_fname);
    
    // Getters

    template <typename T>
    T GetValue(std::string key);

    template <typename T>
    std::vector<T> GetVector(std::string key);

    template <typename T, int N>
    std::array<T, N> GetArray(std::string key);

    inline std::map<std::string, double> GetDictionary(std::string key);

    // Methods

    // void Validate(void);
    
    // I/O

    // friend std::ostream& operator<<(std::ostream& os, const GeneralParameters& _gp);  
};


//
// Constructor
//

inline GeneralParameters::GeneralParameters(const char *_fname) : fname(_fname), table(cpptoml::parse_file(_fname)) {
}

template <typename T>
inline T GeneralParameters::GetValue(std::string key) {
    return table->get_qualified_as<T>(key).value_or(T());
}

template <>
inline int GeneralParameters::GetValue<int>(std::string key) {
    return static_cast<int>(GetValue<double>(key));
}

template <typename T>
inline std::vector<T> GeneralParameters::GetVector(std::string key) {
    return table->get_qualified_array_of<T>(key).value_or(std::vector<T>());
}

template <>
inline std::vector<int> GeneralParameters::GetVector<int>(std::string key) {
    std::vector<int64_t> v = GetVector<int64_t>(key);
    std::vector<int> vout(v.begin(), v.end());
    return vout;
}

// template <typename T>
inline std::map<std::string, double> GeneralParameters::GetDictionary(std::string key) {
    std::map<std::string, double> vout;

    std::shared_ptr<cpptoml::table> inner = table->get_table_qualified(key);
    
    for (auto v: *inner) {
        std::string str = v.first;
        double val = inner->get_qualified_as<double>(str).value_or(0);
        vout.insert(std::make_pair(str, val));
    }

    return vout;
}

//
// I/O
//

// std::ostream& operator<<(std::ostream& os, const GeneralParameters _gp) {  
//     os << "";  
//     return os;  
// }  


#endif  // IO_PARSER_H_
