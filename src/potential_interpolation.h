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
 
class Coordinate_2d;
 
namespace Cartesian_2d
{

class Potential_2d
{
    public:
        virtual double operator() (double x, double y) const  = 0;
};
typedef std::shared_ptr<Potential_2d> Sptr_Potential_2d;

/* Interpolation Potential */

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

//namespace Cylindrical
//{
///* base Potential class */
//class Potential:
//{
//    public:
//        virtual double operator()(double rho, double z, double phi=0) const = 0;
//};

///* potential from Cartesian_2d to Cylindrical */
//class Potential_rotate: Potential
//{
//    public:        
//        Potential_rotate(const Cartesian_2d::Sptr_Potential_2d& pphi)
//            : _pphi(pphi) {}
//        double operator()(double rho, double z, double phi=0) const
//        {
//            return _pphi(rho, z);
//        }
//    private:
//        Cartesian_2d::Sptr_Potential_2d _pphi;    
//};

//} // namespace Cylindrical

namespace Cartesian_3d
{

/* potential from Cylindrical to Cartesian_3d */
class Potential_rotate: public Potential_3d
{
    public:
        Potential_rotate(std::shared_ptr<Cartesian_2d::Potential_2d> pphi)
            : _pphi(pphi) {}   
        double operator() (double x, double y, double z) const
        {
            return (*_pphi)(gsl_hypot(x, y), z);
        }
    private:    
        const std::shared_ptr<Cartesian_2d::Potential_2d> _pphi;
};


} // namespace Cartesian_3d






#endif // _POTENTIAL_INTERPOLATION_H_
