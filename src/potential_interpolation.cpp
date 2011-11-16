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
using namespace CGAL;

typedef Data_access<Point_value_map> Data_access_value;
typedef Data_access<Point_vector_map> Data_access_grad;

namespace Cartesian_2d
{

/* Interpolation Potential */
Potential_interpolation_2d::Potential_interpolation_2d(
    const Coordinate_2d& coords, const vector<double>& values)
{
    size_t n = coords.size();
    for (size_t i = 0; i < n; ++i)
    {
        Point p(coords.x(i), coords.y(i));
        _dt.insert(p);
        _value_map.insert(make_pair(p, values[i]));
    }
    sibson_gradient_fitting_nn_2(
        _dt, inserter(_grad_map, _grad_map.begin()),
        Data_access_value(_value_map),
        GradTraits());
}

double Potential_interpolation_2d::operator()(double rho, double z) const
{
    Point p(rho, z);
    vector<pair<Point, double> > ncoords;
    double norm = natural_neighbor_coordinates_2(
        _dt, p, back_inserter(ncoords)).second;
    pair<double, bool> res = quadratic_interpolation(
        ncoords.begin(), ncoords.end(), norm, p,
        Data_access_value(_value_map),
        Data_access_grad(_grad_map),
        GradTraits());
    if (res.second)
    {
        return res.first;
    }
    else
    {
        return linear_interpolation(
            ncoords.begin(), ncoords.end(), norm,
            Data_access_value(_value_map));
    }
}

} // namespace Cartesian_2d
