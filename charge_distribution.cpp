#include <config.h>

#include <cmath>
#include <cassert>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_linalg.h>


static const double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;

gsl_vector* Charge_distri(double phi, double R, const gsl_matrix* inner_dist)
/* input: potential, radius, inner_distance */
{
    size_t n = inner_dist->size1;
    gsl_matrix* A = gsl_matrix_alloc(n, n);
    for (size_t i = 0; i < n; ++i)
    {   for (size_t j = 0; j < i; ++j)
        {
            double rr = gsl_matrix_get(inner_dist, i, j);
            rr = 1 / rr;
            gsl_matrix_set(A, i, j, rr);
            gsl_matrix_set(A, j, i, rr);
        }
    }
    
    gsl_vector_view diag = gsl_matrix_diagonal(A);
    gsl_vector_set_all(&diag.vector, 1/R);
    
    assert(gsl_matrix_max(A) == 1/R);
    
    gsl_matrix_scale(A, 1/(4*M_PI*epsilon_0));
    
    gsl_permutation* p = gsl_permutation_alloc(n);    
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);
    
    gsl_vector* b = gsl_vector_alloc(n);
    gsl_vector_set_all(b, phi);
    
    gsl_vector* x = gsl_vector_alloc(n);
    gsl_linalg_LU_solve(A, p, b, x);
    
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_permutation_free(p);
    
    return x;
}
