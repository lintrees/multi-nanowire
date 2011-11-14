#include <config.h>

#include <vector>
#include <string>
#include <sstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#include "coordinate.h"


Coordinate_2d::Coordinate_2d(std::istream& is)
{
    std::string str;
    std::stringstream ss;
    std::vector<double> vx, vy;
    double x, y;
    while(getline(is, str))
    {
        ss.clear();
        ss.str(str);
        ss >> x >> y;
        vx.push_back(x);
        vy.push_back(y);
    }
    size_t n = vx.size();
    gsl_vector_view vvx = gsl_vector_view_array(vx.data(), n);
    gsl_vector_view vvy = gsl_vector_view_array(vy.data(), n);
    _m = gsl_matrix_alloc(n, 2);
    gsl_matrix_set_col(_m, 0, &vvx.vector);
    gsl_matrix_set_col(_m, 1, &vvy.vector);
}
Coordinate_2d::~Coordinate_2d()
{
    gsl_matrix_free(_m);
}

Inner_distance::Inner_distance(const Coordinate_2d& coord)
{
    size_t n = coord.size();
    _m = gsl_matrix_alloc(n, n);    
    gsl_matrix_const_view
        coord_view = gsl_matrix_const_view_array(coord.data(), n, 2);
    gsl_vector_const_view xv = gsl_matrix_const_column(&coord_view.matrix, 0);
    gsl_vector_const_view yv = gsl_matrix_const_column(&coord_view.matrix, 1);
    for (size_t i=0, j=1; j < n; ++i, ++j)
    {
        gsl_vector_view xv1 = gsl_matrix_subcolumn(_m, i, j, n-j);
        gsl_vector_const_view xv2 = gsl_vector_const_subvector(&xv.vector, j, n-j);
        double mx2 = -gsl_vector_get(&xv.vector, i);
        gsl_vector_memcpy(&xv1.vector, &xv2.vector);
        gsl_vector_add_constant(&xv1.vector, mx2);

        gsl_vector_view yv1 = gsl_matrix_subrow(_m, i, j, n-j);
        gsl_vector_const_view yv2 = gsl_vector_const_subvector(&yv.vector, j, n-j);
        double my2 = -gsl_vector_get(&yv.vector, i);
        gsl_vector_memcpy(&yv1.vector, &yv2.vector);
        gsl_vector_add_constant(&yv1.vector, my2);
    }
    for (size_t i = 0; i < n; ++i)
    {   for (size_t j = i+1; j < n; ++j)
        {
            double dy = gsl_matrix_get(_m, i, j);
            double dx = gsl_matrix_get(_m, j, i);
            double dist = gsl_hypot(dx, dy);
            gsl_matrix_set(_m, i, j, dist);
            gsl_matrix_set(_m, j, i, dist);
        }
    }
}

Inner_distance::~Inner_distance()
{
    gsl_matrix_free(_m);
}
