#include <config.h>

#include <iostream>
#include <vector>
#include <cassert>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_multifit_nlin.h>

#include <cstdio>
#include <cstdlib>

static constexpr double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
static const double coeff = 1/ (epsilon_0 * 4*M_PI);

static std::vector<double> X, Y;
static constexpr double epsabs = 1e10;

static constexpr size_t p = 3;
static size_t n;


static int f_f(const gsl_vector* x, void* params, gsl_vector* f)
{
    double R = gsl_vector_get(x, 0);
    double Q = gsl_vector_get(x, 1);
    double phi0 = gsl_vector_get(x, 1);
    for (size_t i = 0; i < n; ++i)
    {
        double xi = X[i], yi = Y[i];     
        double d = xi+R;   
        gsl_vector_set(f, i, coeff* Q/d -phi0 -yi );
    }
    return 0;
}

static int f_df(const gsl_vector* x, void* params, gsl_matrix* J)
{
    double R = gsl_vector_get(x, 0);
    double Q = gsl_vector_get(x, 1);
    double phi0 = gsl_vector_get(x, 1);
    for (size_t i = 0; i < n; ++i)
    {
        double xi = X[i], yi = Y[i];
        double d = xi+R;
        gsl_matrix_set(J, i, 0, coeff* -Q/(d*d) );
        gsl_matrix_set(J, i, 1, coeff* 1/d );
        gsl_matrix_set(J, i, 2, -1 );
    }
    return 0;
}

static int f_fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
{
    double R = gsl_vector_get(x, 0);
    double Q = gsl_vector_get(x, 1);
    double phi0 = gsl_vector_get(x, 1);
    for (size_t i = 0; i < n; ++i)
    {
        double xi = X[i], yi = Y[i];
        double d = xi+R;   
        gsl_vector_set(f, i, coeff* Q/d -phi0 -yi );
        gsl_matrix_set(J, i, 0, coeff* -Q/(d*d) );
        gsl_matrix_set(J, i, 1, coeff* 1/d );
        gsl_matrix_set(J, i, 2, -1 );
    }
    return 0;
}

int main(int argc, char** argv)
{
//    assert(argc==2);
//    phi0 = atof(argv[1]);
    double xi, yi;
    while (std::cin >> xi >>yi)
    {
        X.push_back(xi);
        Y.push_back(yi);        
    }
    n = X.size();

    gsl_multifit_fdfsolver* 
        s = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, n, p);
    gsl_multifit_function_fdf f;
    f.f = &f_f;
    f.df = &f_df;
    f.fdf = &f_fdf;
    f.n = n;
    f.p = p;
    
    double x_init[p] = {1e-8, -1e-17, -10};
    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);

    gsl_vector* g = gsl_vector_alloc(p);
    int status;
    do
    {
        status = gsl_multifit_fdfsolver_iterate(s);
        assert(!status); // FIXME 
        printf(".");
        if (status)        
        {   printf("\n");break; }
        gsl_multifit_gradient(s->J, s->f, g);
        status = gsl_multifit_test_gradient(g, epsabs);
    } while (status == GSL_CONTINUE);
    
    double R = gsl_vector_get(s->x, 0);
    printf("%.12g\n", R);
    double Q = gsl_vector_get(s->x, 1);
    printf("%.12g\n", Q);
    double phi0 = gsl_vector_get(s->x, 2);
    printf("%.12g\n", phi0);
}


