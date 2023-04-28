#ifndef DATA_HPP
#define DATA_HPP


#include <string>
#include <unordered_map>

#include <boost/any.hpp>
#include <boost/unordered_map.hpp>


class Data
{

private:

    boost::unordered_map<std::string, boost::any> _data;


public:

    Data(void);
    Data(std::string, boost::any value);

    void SetValue(std::string key, boost::any value);

    template <typename T>
    void IncrementValue(std::string key, T value);

    template <typename T>
    T GetValue(std::string key);

};


inline Data::Data(void)
{}

inline Data::Data(std::string key, boost::any value)
: _data({{key, value}})
{}

inline void Data::SetValue(std::string key, boost::any value)
{
    _data[key] = value;
}

template <typename T>
inline void Data::IncrementValue(std::string key, T value)
{
    _data[key] = boost::any_cast<T>(_data[key]) + value;
}

template <typename T>
inline T Data::GetValue(std::string key)
{
    return boost::any_cast<T>(_data[key]);
}


#endif // DATA_HPP