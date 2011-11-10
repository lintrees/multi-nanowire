#ifndef _POTENTIAL_INTERPOLATION_H_
#define _POTENTIAL_INTERPOLATION_H_

#include "potential.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation;

typedef CGAL::Interpolation_gradient_fitting_traits_2<K> GradTraits;
typedef K::Point_2 Point;
typedef std::map<Point, K::FT, K::Less_xy_2> Point_value_map;
typedef std::map<Point, K::Vector_2 , K::Less_xy_2> Point_vector_map;


/* 
 * 2d Potential 
 */
 
namespace Cartesian_2d
{

class Potential_2d
{
    public:
        virtual double operator() (double x, double y) const  = 0;
};

/* Interpolation Potential */
class Coordinate_2d;
class Potential_interpolation_2d: Potential_2d
{
    public:
        Potential_interpolation_2d(const Coordinate_2d& , const std::vector<double>&);
        double operator()(double x, double y) const;
    private:
        Delaunay_triangulation _dt;
        Point_value_map _value_map;
        Point_vector_map _grad_map;
};

} // namespace Cartesian_2d

namespace Cylindrical
{

class Potential:
{
    public:
        virtual double operator()(double rho, double z, double phi=0) const  = 0;
};

class Potential_section:Potential
{
    
};


} // namespace Cylindrical


#endif // _POTENTIAL_INTERPOLATION_H_
