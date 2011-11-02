#include <config.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#include <cstdio>

const char* fn_array_coord = "array_coord.dat";

struct Nanowire
{
    double x, y;
};


void Readin_coord(std::vector<Nanowire>& nanowires)
{
    std::ifstream ifs(fn_array_coord);
    assert(ifs);
    
    Nanowire nanowire;
    double x, y;
    while (ifs >> x >> y)
    {
        nanowire = Nanowire();
        nanowire.x = x;
        nanowire.y = y;
        nanowires.push_back(nanowire);
    }
}

gsl_matrix* Coord_matrix(const std::vector<Nanowire>& nanowires)
{
    size_t n = nanowires.size();
    gsl_matrix* coord_matrix = gsl_matrix_alloc(n, 2);
    size_t i = 0;
    for (std::vector<Nanowire>::const_iterator ptr = nanowires.begin(); ptr != nanowires.end(); ++ptr, ++i)
    {
        gsl_matrix_set(coord_matrix, i, 0, ptr->x);
        gsl_matrix_set(coord_matrix, i, 1, ptr->y);
    }
    return coord_matrix;
}

gsl_matrix* Inner_distance(const gsl_matrix* coord)
{
    size_t n = coord->size1;
    gsl_vector_const_view xv = gsl_matrix_const_column(coord, 0);
    gsl_vector_const_view yv = gsl_matrix_const_column(coord, 1);
    gsl_matrix* m = gsl_matrix_calloc(n, n); // inner distance matrix

    for (size_t i=0, j=1; j < n; ++i, ++j)
    {
        gsl_vector_view xv1 = gsl_matrix_subcolumn(m, i, j, n-j);
        gsl_vector_const_view xv2 = gsl_vector_const_subvector(&xv.vector, j, n-j);
        double mx2 = -gsl_vector_get(&xv.vector, i);
        gsl_vector_memcpy(&xv1.vector, &xv2.vector);
        gsl_vector_add_constant(&xv1.vector, mx2);

        gsl_vector_view yv1 = gsl_matrix_subrow(m, i, j, n-j);
        gsl_vector_const_view yv2 = gsl_vector_const_subvector(&yv.vector, j, n-j);
        double my2 = -gsl_vector_get(&yv.vector, i);
        gsl_vector_memcpy(&yv1.vector, &yv2.vector);
        gsl_vector_add_constant(&yv1.vector, my2);
    }

    for (size_t i = 0; i < n; ++i)
    {   for (size_t j = i+1; j < n; ++j)
        {
            double dy = gsl_matrix_get(m, i, j);
            double dx = gsl_matrix_get(m, j, i);
            double dist = gsl_hypot(dx, dy);
            gsl_matrix_set(m, i, j, dist);
            gsl_matrix_set(m, j, i, dist);
        }
    }

    return m;
}

