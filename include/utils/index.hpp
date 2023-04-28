#ifndef UTILS_INDEX_H_
#define UTILS_INDEX_H_


#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>

#include "xtensor/xarray.hpp"


//
// Convert linear index to subscripts
//

static inline std::vector<int> ind2sub(const int i, const std::vector<int> n) {
    int l = i;
    std::vector<int> loc(n.size(), 0);
    std::vector<int> v = n;

    std::partial_sum(v.begin(), v.end(), v.begin(), std::multiplies<int>());
    std::rotate(v.rbegin(), v.rbegin() + 1, v.rend());
    v[0] = 1;

    for (int ii = loc.size()-1; ii >= 0; --ii) {
        std::div_t x = std::div(l, v[ii]);
        loc[ii] = x.quot;
        l = x.rem;
    }

    return loc;
}

//
// Convert subscripts to linear index
//

inline int sub2ind(const std::vector<int> i, const std::vector<int> n) {
    int loc = 0;
    std::vector<int> v = n;
    
    std::partial_sum(v.begin(), v.end(), v.begin(), std::multiplies<int>());
    std::rotate(v.rbegin(), v.rbegin() + 1, v.rend());
    v[0] = 1;
    
    return std::inner_product(v.begin(), v.end(), i.begin(), 0);
}

inline int sub2ind(xt::xarray<int> i, xt::xarray<int> n)
{
    return sub2ind(
        std::vector<int>(i.begin(), i.end()), 
        std::vector<int>(n.begin(), n.end())
    );
}

inline int sub2ind(std::vector<int> i, xt::xarray<int> n)
{
    return sub2ind(
        std::vector<int>(i.begin(), i.end()), 
        std::vector<int>(n.begin(), n.end())
    );
}

//
// Create neighbour vectors
//
/*
    This function generates a face neighbouring vector starting at the
    negative and ending at the positive face. 

    - 1d: (-1)      , (+1)
    - 2d: (-1, 0)   , (+1, 0),  , (0, -1)   , ...
    - 3d: (-1, 0, 0), (+1, 0, 0), (0, -1, 0), ...
*/

inline xt::xarray<int> direction_vector(void) 
{
    return xt::xarray<int> {
        {-1,  0,  0},
        {+1,  0,  0},
        { 0, -1,  0},
        { 0, +1,  0},
        { 0,  0, -1},
        { 0,  0, +1}
    };
}

inline std::vector<int> face_vectors(int dimension) {

    std::vector<int> neighbours = 
        std::vector<int>(2*std::pow(dimension, 2), 0);

    for (int j = 0; j < neighbours.size(); ++j) {
        std::vector<int> face = std::vector<int>(dimension, 0);
        face[0]               = 1;

        int irot = j % dimension;
        int ipow = j / dimension;

        std::rotate(face.rbegin(), face.rbegin() + static_cast<int>(ipow/2), face.rend());
        neighbours[j] = face[irot] * (std::pow(-1, ipow + 1));
    }

    return neighbours;
}

inline std::vector<int> face_neighbours_index(int _dimension, int _index, std::vector<int> _shape) {
    
    std::vector<int> neighbours(2*_dimension);
    std::vector<int> face_vector = face_vectors(_dimension);
    std::vector<int> subind      = ind2sub(_index, _shape);

    for (int i=0; i<neighbours.size(); ++i) {
        std::vector<int> temp = 
            std::vector<int>(&face_vector[3*i], &face_vector[3*(i+1)]);

        std::transform(
            subind.begin(),
            subind.end(),
            temp.begin(),
            temp.begin(),
            std::plus<double>()
        );

        neighbours[i] = sub2ind(temp, _shape);
    }

    return neighbours;
}

inline std::vector<int> face_neighbours_no_boundary_index(int _dimension, int _index, std::vector<int> _shape) {
        
    std::vector<int> neighbours(2*_dimension, -1);
    std::vector<int> face_vector = face_vectors(_dimension);
    std::vector<int> subind      = ind2sub(_index, _shape);

    for (int i=0; i<neighbours.size(); ++i) {
        std::vector<int> temp = 
            std::vector<int>(&face_vector[_dimension*i], &face_vector[_dimension*(i+1)]);

        std::transform(
            temp.begin(),
            temp.end(),
            subind.begin(),
            temp.begin(),
            std::plus<double>()
        );

        bool is_equal = false;
        for (int j=0; j<temp.size(); ++j) {
            if (temp[j] == _shape[j] || temp[j] == -1) {
                is_equal = true;
            }
        }

        if (is_equal) continue;

        neighbours[i] = sub2ind(temp, _shape);
    }

    std::vector<int>::iterator pend = 
        std::remove_if(neighbours.begin(), neighbours.end(), [](int i) {return i==-1;});

    return std::vector<int>(neighbours.begin(), pend);
}

inline bool is_boundary(std::vector<int> _sub_index, std::vector<int> _shape) {

    for (int i=0; i<_sub_index.size(); ++i) {
        if (_sub_index[i] >= _shape[i] - 2) return true;
        if (_sub_index[i] <= 1            ) return true;
    }

    return false;
}

inline bool is_boundary(int _index, std::vector<int> _shape) {
    return is_boundary(ind2sub(_index, _shape), _shape);
}

inline std::vector<int> mirror_subindex(std::vector<int> _sub_index, std::vector<int> _shape, int _lag=1) {
    std::vector<int> out(_sub_index);

    for (int i=0; i<_shape.size(); ++i) {
        if (_sub_index[i] == 0          ) out[i] = _shape[i]-1-_lag;
        if (_sub_index[i] == _shape[i]-1) out[i] = _lag;
    }

    return out;
}

inline int mirror_index(int _index, std::vector<int> _shape, int _lag=1) {
    return sub2ind(mirror_subindex(ind2sub(_index, _shape), _shape, _lag), _shape);
}


#endif  // UTILS_INDEX_H_
