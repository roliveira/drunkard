#ifndef UTILS_INTERPOLATION_H_
#define UTILS_INTERPOLATION_H_

#include <vector>
#include <array>
/*
#include "linterp.h"

#include "utils/root.hpp"


inline std::unique_ptr<NDInterpolator_1_S> interpolate(std::vector<double> xvalues, std::vector<double> yvalues) {
    std::vector< std::vector<double>::iterator> grid_iter_list;
    grid_iter_list.push_back(xvalues.begin());

    std::array<int, 1> grid_size = { static_cast<int>(xvalues.size()) };
    int num_elements = grid_size[0];

    std::unique_ptr<NDInterpolator_1_S> interp_simplex(
        new NDInterpolator_1_S(
            grid_iter_list.begin(), 
            grid_size.begin(), 
            yvalues.data(), 
            yvalues.data() + num_elements
        )
    );

    return std::move(interp_simplex);
}
*/

#endif  // UTILS_INTERPOLATION_H_
