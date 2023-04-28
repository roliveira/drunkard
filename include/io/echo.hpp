#ifndef IO_ECHO_H_
#define IO_ECHO_H_


#include <iostream>
#include <fstream>
#include <string>
#include <sstream> 
#include <vector>
#include <iomanip>

#include "config.hpp"


inline std::string version(void) {
    return std::to_string(DRUNKARD_VERSION_MAJOR) + "." 
         + std::to_string(DRUNKARD_VERSION_MINOR) + "." 
         + std::to_string(DRUNKARD_VERSION_PATCH);
}

inline std::string banner(void) {
    std::string NAME = "drunkard";
    std::string VERSION = version();
    std::string BANNER = "";

    int num_spaces = NAME.size() + VERSION.size();

    BANNER += "\u250C"; for(int i = 0; i < NUM_CHAR_MESSAGE-2; ++i) BANNER += "\u2500"; BANNER += "\u2510\n";
    BANNER += "\u2502" + std::string(NUM_CHAR_MESSAGE-2, ' ') + "\u2502\n";
    
    BANNER += "\u2502" + std::string((NUM_CHAR_MESSAGE-num_spaces)/2-1, ' ');
    BANNER += NAME + " " + VERSION;
    BANNER += std::string((NUM_CHAR_MESSAGE-num_spaces)/2-1, ' ') + "\u2502\n";

    BANNER += "\u2502" + std::string(NUM_CHAR_MESSAGE-2, ' ') + "\u2502\n";
    BANNER += "\u2514"; for(int i = 0; i < NUM_CHAR_MESSAGE-2; ++i) BANNER += "\u2500"; BANNER += "\u2518\n";

    return BANNER;
}

inline std::string separator(void) {
    std::string SEPARATOR = "\n";
    for(int i = 0; i < NUM_CHAR_MESSAGE; ++i) SEPARATOR += "\u2500";
    SEPARATOR += "\n";
    return SEPARATOR;
}

inline std::string section(std::string title) {
    std::string out = title + " ";
    int len = out.size();
    for(int i = 0; i < NUM_CHAR_MESSAGE-len; ++i) out += "\u2500";
    out = "\n\n" + out + "\n\n";
    return out;
}

inline std::string section_close(std::string title, double time) {
    std::stringstream out;
    std::stringstream val;
    int len;
    int time_trunc;
    int hh, mm, ss;

    time_trunc = static_cast<int>(time);
    hh         = time_trunc / 3600;
	time_trunc = time_trunc % 3600;
	mm         = time_trunc / 60;
	time_trunc = time_trunc % 60;
	ss         = time_trunc;
 
    out << title + ' ';

    val << " " << std::setw(2) << std::setfill('0') << hh;
    val << ":" << std::setw(2) << std::setfill('0') << mm;
    val << ":" << std::setw(2) << std::setfill('0') << ss;
    
    len = out.tellp() + val.tellp();
    for(int i = 0; i < NUM_CHAR_MESSAGE-len; ++i) out << "\u2500";
    out << val.rdbuf();

    return out.str();
}

template <typename T>
inline std::string message(std::string info, std::vector<T> value) {
    std::stringstream out;
    std::stringstream val;
    int len;

    out << info << ' ';
    
    val << std::showpos << std::scientific << *value.cbegin();
    for (typename std::vector<T>::const_iterator it = value.cbegin()+1; it != value.cend(); ++it) {
        val << std::showpos << std::scientific << ", " << *it;
    }

    len = out.tellp() + val.tellp();

    if (len < NUM_CHAR_MESSAGE) {
        out << std::string(NUM_CHAR_MESSAGE - len, ' ');
    }

    out << val.rdbuf() << std::endl;

    return out.str();
}

template <typename T>
inline std::string message(std::string info, T value) {
    std::stringstream out;
    std::stringstream val;
    int len;

    out << info << ' ';
    val << std::showpos << std::scientific << value;
    len = out.tellp() + val.tellp();

    if (len < NUM_CHAR_MESSAGE) {
        out << std::string(NUM_CHAR_MESSAGE - len, ' ');
    }

    out << val.rdbuf() << std::endl;

    return out.str();
}


#endif  // IO_ECHO_H_
