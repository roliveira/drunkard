#ifndef UTILS_OPERATIONS_H
#define UTILS_OPERATIONS_H


#include <numeric>
#include <cmath>
#include <vector>

#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xio.hpp"


inline double VectorMagnitude(std::vector<double> source_position, std::vector<double> target_position) {
    double magnitude = 0;
    
    for (int i = 0; i < source_position.size(); ++i) 
        magnitude += std::pow(target_position.at(i) - source_position.at(i), 2);

    return std::sqrt(magnitude);
}

inline double VectorMagnitude(xt::xarray<double> source_position, xt::xarray<double> target_position)
{
    return xt::sqrt(xt::sum(xt::pow(target_position-source_position, 2)))(0);
}

/*inline std::array<double, 3> VectorDirection(std::array<double, 3> source_position, std::array<double, 3> target_position) {
    double magnitude = VectorMagnitude(source_position, target_position);
    
    std::array<double, 3> direction;
    for (int i = 0; i < source_position.size(); ++i) direction[i] = (target_position[i] - source_position[i]) / magnitude;

    return direction;
}

inline std::array<double, 3> VectorReflection(std::array<double, 3> incident_direction, std::array<double, 3> normal_direction) {
    double prod = std::inner_product(incident_direction.begin(), incident_direction.begin(), normal_direction.begin(), 0.0);

    std::array<double, 3> reflection;
    for (int i = 0; i < incident_direction.size(); ++i) reflection[i] = 2 * prod * normal_direction[i] - incident_direction[i];

    return reflection;
}

inline std::array<double, 3> PointOfIntersection(std::array<double, 3> point, std::array<double, 3> vector_direction, std::array<double, 3> face_normal, std::array<double, 3> face_centre) {
    double constant = std::inner_product(face_normal.begin(), face_normal.end(), face_centre.begin(), 0.0);
    double t = (constant - std::inner_product(face_normal.begin(), face_normal.end(), point.begin(), 0.0))
        / std::inner_product(face_normal.begin(), face_normal.end(), vector_direction.begin(), 0.0);

    std::array<double, 3> interception;
    for (int i = 0; i < interception.size(); ++i) interception[i] = vector_direction[i] * t + point[i];

    return interception;
}*/

#endif  // UTILS_OPERATIONS_H
