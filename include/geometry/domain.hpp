
#ifndef GEOMETRY_DOMAIN_H_
#define GEOMETRY_DOMAIN_H_

#include <vector>
#include <random>

#include "utils/index.hpp"
#include "geometry/geometry.hpp"


class DomainSystem : public GeometrySystem {
private:

    int _label_background;

    std::vector<int> _label;
    std::vector<int> _label_unique;

public:

    // Constructors

    DomainSystem(void);
    DomainSystem(std::vector<int> shape);
    DomainSystem(std::vector<int> shape, int label_background);

    // Getters

    const std::vector<int> GetLabel(void);
    const int GetLabel(int index);

    const std::vector<int> GetLabelUnique(void);

    // Setters

    void SetLabel(std::vector<int> values);
    void SetLabel(int value, int index);
    
    // void SetLabelBackground(int label_background);
    // void SetLabelUnique(std::vector<int> values);

    // Methods

    void CreateRectangle(int label, std::vector<int> bounding_coordinates);
    void CreateCircle(int label, std::vector<int> center, int radius);
    void CreateRandom(std::vector<int> labels);
};

// 
// Constructors
// 

inline DomainSystem::DomainSystem(void)
    : GeometrySystem() 
{}

inline DomainSystem::DomainSystem(std::vector<int> shape)
    : GeometrySystem(shape) {
    _label = std::vector<int>(this->GetLinearShape(), 0);
}

inline DomainSystem::DomainSystem(std::vector<int> shape, int label_background)
    : GeometrySystem(shape), _label_background(label_background) {
    _label = std::vector<int>(this->GetLinearShape(), _label_background);
}

// 
// Getters
// 

inline const std::vector<int> DomainSystem::GetLabel(void) {
    return _label;
}

inline const int DomainSystem::GetLabel(int index) {
    return _label.at(index);
}

inline const std::vector<int> DomainSystem::GetLabelUnique(void) {
    return _label_unique;
}

// 
// Setters
// 

inline void DomainSystem::SetLabel(std::vector<int> values) {
    _label = values;
}

inline void DomainSystem::SetLabel(int value, int index) {
    _label[index] = value;
}

// 
// Methods
// 

// inline std::vector<int> DomainSystem::SetLabelUnique(std::vector<int> values) {}

/**
 * @brief Creates a rectangle using the bounding_coordinates as a limit and assings
 * label to it. 
 * 
 * @param label                : Label value inside the rectangle
 * @param bounding_coordinates : Coordinates as [x0, y0, x1, y1] index values. 
 */
inline void DomainSystem::CreateRectangle(int label, std::vector<int> bounding_coordinates) {
    for (int i=0; i<this->GetLinearShape(); ++i) {
        std::vector<int> subind = ind2sub(i, this->GetShape());

        if ((subind[0] >= bounding_coordinates[0]) 
         && (subind[1] >= bounding_coordinates[1])
         && (subind[0] <= bounding_coordinates[2])
         && (subind[1] <= bounding_coordinates[3]))
            SetLabel(label, i);
    }
}

/**
 * @brief Creates a circle defined by its center and radius and assigns a label to it.
 * 
 * @param label  : Label value inside the circle
 * @param center : Center of the circle as [xc, yc] index values
 * @param radius : Radius of the circle as an index value.
 */
inline void DomainSystem::CreateCircle(int label, std::vector<int> center, int radius) {
    for (int i=0; i<this->GetLinearShape(); ++i) {
        std::vector<int> subind = ind2sub(i, this->GetShape());

        if (std::pow(subind[0]-center[0], 2) + std::pow(subind[1]-center[1], 2) <= std::pow(radius, 2))
            SetLabel(label, i);
    }
}

/**
 * @brief Randomly assigns a label from labels to the system.
 * 
 * @param labels : The list of labels to choose from [0, 1, ..., N]
 */
inline void DomainSystem::CreateRandom(std::vector<int> labels) {
    std::mt19937 generator;
    std::uniform_int_distribution<int> dist(0, labels.size()-1);

    for (int i=0; i<this->GetLinearShape(); ++i) {
        SetLabel(labels.at(dist(generator)), i);
    }
}

#endif  // GEOMETRY_DOMAIN_H_
