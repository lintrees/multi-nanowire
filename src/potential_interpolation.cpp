#include <config.h>

#include <cmath>
#include <cassert>
#include <iostream>

#include <gsl/gsl_const.h>

#include "potential_interpolation.h"
#include "coordinate.h"

#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/sibson_gradient_fitting.h>
#include <CGAL/interpolation_functions.h>

using namespace std;

namespace Cartesian_2d
{

/* Interpolation Potential */
Potential_interpolation_2d::Potential_interpolation_2d(
    const Coordinate_2d& coords, const std::vector<double>& values)
{
    size_t n = coords.size();
    for (size_t i = 0; i < n; ++i)
    {
        Point p(coords.x(i), coords.y(i));
        _dt.insert(p);
        _value_map.insert(std::make_pair(p, values[i]));
    }    
    CGAL::sibson_gradient_fitting_nn_2(
        _dt, std::inserter(_grad_map, _grad_map.begin()),
        CGAL::Data_access<Point_value_map>(_value_map),
        GradTraits());
}

double Potential_interpolation_2d::operator()(double rho, double z) const
{
    Point p(rho, z);
    std::vector<std::pair<Point, double> > ncoords;
    double norm = CGAL::natural_neighbor_coordinates_2(
        _dt, p, std::back_inserter(ncoords)).second;        
    std::pair<double, bool> res = CGAL::quadratic_interpolation(
        ncoords.begin(), ncoords.end(), norm, p,
        CGAL::Data_access<Point_value_map>(_value_map),
        CGAL::Data_access<Point_vector_map>(_grad_map),
        GradTraits());        
    if (!res.second)
    {
        return CGAL::linear_interpolation(
            ncoords.begin(), ncoords.end(), norm,
            CGAL::Data_access<Point_value_map>(_value_map));
    }
    assert(res.second);
    return res.first;
}

} // namespace Cartesian_2d
