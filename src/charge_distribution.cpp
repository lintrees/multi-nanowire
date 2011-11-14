#include <config.h>

#include <cmath>
#include <cassert>
#include <vector>
#include <utility>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_linalg.h>

#include "coordinate.h"

using namespace std;

static constexpr double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;

vector<double> Charge_distri(const vector<double>& vphi, double R, const Inner_distance& inner_dist)
/* input: potential, radius, inner_distance */
{
    size_t n = inner_dist.size();
    gsl_matrix* A = gsl_matrix_alloc(n, n);
    double R2 = 2*R;
    for (size_t i = 0; i < n; ++i)
    {   for (size_t j = 0; j < i; ++j)
        {
            double rr = inner_dist(i, j);
            assert(rr > R2);
            rr = 1 / rr;
            gsl_matrix_set(A, i, j, rr);
            gsl_matrix_set(A, j, i, rr);
        }
    }
    gsl_vector_view diag = gsl_matrix_diagonal(A);
    gsl_vector_set_all(&diag.vector, 1/R);
    gsl_matrix_scale(A, 1/(4*M_PI*epsilon_0));

    gsl_permutation* p = gsl_permutation_alloc(n);
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);
    
    gsl_vector_const_view  vb = 
        gsl_vector_const_view_array(vphi.data(), n);
    gsl_vector* x = gsl_vector_alloc(n);
    gsl_linalg_LU_solve(A, p, &vb.vector, x);

    std::vector<double> res(x->data, x->data+n);

    gsl_matrix_free(A);
    gsl_vector_free(x);
    gsl_permutation_free(p);

    return std::move(res);
}
