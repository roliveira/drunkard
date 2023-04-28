#ifndef PARTICLE_CELL_HPP_
#define PARTICLE_CELL_HPP_


#include <list>

#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xio.hpp"

#include "utils/index.hpp"
#include "particle.hpp"


class Cell
{

public:

    xt::xarray<int> _num_grid;  //!< Number of grid blocks in a cell
    xt::xarray<int> _shape_grid;
    xt::xarray<int> _shape_cell;

    xt::xarray<std::list<Particle*>> _parts;
    
    Cell(void);
    Cell(std::vector<int> num_grid, std::vector<int> shape_grid);

    void push_back(Particle* p, int ci, int cj, int ck);
    void push_front(Particle* p, int ci, int cj, int ck);

    std::vector<int> Grid2Cell(std::vector<int> grid_index);
    int Sub2Ind(std::vector<int> subind);
};


Cell::Cell(void)
{}

Cell::Cell(std::vector<int> num_grid, std::vector<int> shape_grid)
{
    this->_num_grid = xt::adapt(num_grid, {3});
    this->_shape_grid = xt::adapt(shape_grid, {3});
    this->_shape_cell = this->_shape_grid / this->_num_grid;

    this->_parts = xt::empty<std::list<Particle*>>(this->_shape_cell);
    this->_parts.fill(std::list<Particle*>());
}

void Cell::push_back(Particle* p, int ci, int cj, int ck)
{
    this->_parts(ci, cj, ck).push_back(p);
}

void Cell::push_front(Particle* p, int ci, int cj, int ck)
{
    this->_parts(ci, cj, ck).push_front(p);
}

std::vector<int> Cell::Grid2Cell(std::vector<int> grid_index)
{
    std::vector<int> cell_index(3, 0);

    for (int i=0; i<cell_index.size(); ++i)
    {
        cell_index[i] = grid_index[i] / this->_num_grid[i];
    }

    return cell_index;
}

int Cell::Sub2Ind(std::vector<int> subind)
{
    return sub2ind(subind, this->_shape_cell);
}


#endif  // PARTICLE_CELL_HPP_