#include <config.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_multifit_nlin.h>

#include <cstdio>
#include <cstdlib>

static const double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
static const double epsilon_0_4_pi = epsilon_0 * 4*M_PI;

static const char* fnin = "data/barrier.dat";

static std::vector<double> X, Y;
static const double epsabs = 1e-8;

static size_t n;
static const size_t p = 2;
static double phi0;


static int f_f(const gsl_vector* x, void* params, gsl_vector* f)
{
    double R = gsl_vector_get(x, 0);
    double Q = gsl_vector_get(x, 1);
    for (size_t i = 0; i < n; ++i)
    {
        double xi = X[i], yi = Y[i];        
        gsl_vector_set(f, i, Q/(epsilon_0_4_pi*(xi+R))-phi0 -yi );
    }
    return 0;
}

static int f_df(const gsl_vector* x, void* params, gsl_matrix* J)
{
    double R = gsl_vector_get(x, 0);
    double Q = gsl_vector_get(x, 1);
    for (size_t i = 0; i < n; ++i)
    {
        double xi = X[i], yi = Y[i];        
        gsl_matrix_set(J, i, 0, -Q/(epsilon_0_4_pi*(xi+R)*(xi+R)) );
        gsl_matrix_set(J, i, 1, 1/(epsilon_0_4_pi*(xi+R)) );
    }
    return 0;
}

static int f_fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
{
    double R = gsl_vector_get(x, 0);
    double Q = gsl_vector_get(x, 1);
    for (size_t i = 0; i < n; ++i)
    {
        double xi = X[i], yi = Y[i];        
        gsl_vector_set(f, i, Q/(epsilon_0_4_pi*(xi+R))-phi0 -yi );
        gsl_matrix_set(J, i, 0, -Q/(epsilon_0_4_pi*(xi+R)*(xi+R)) );
        gsl_matrix_set(J, i, 1, 1/(epsilon_0_4_pi*(xi+R)) );
    }
    return 0;
}

int main(int argc, char** argv)
{
    assert(argc==2);
    phi0 = atof(argv[1]);
    std::ifstream ifs(fnin);
    assert(ifs);
    double xi, yi;
    while (ifs >> xi >>yi)
    {
        X.push_back(xi);
        Y.push_back(yi);        
    }
    ifs.close();
    n = X.size();

    gsl_multifit_fdfsolver* s = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, n, p);

    gsl_multifit_function_fdf f;
    f.f = &f_f;
    f.df = &f_df;
    f.fdf = &f_fdf;
    f.n = n;
    f.p = p;
    
    double x_init[p] = {1e-8, -1e-16};
    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);

    gsl_vector* g = gsl_vector_alloc(p);
    int status;
    do
    {
        status = gsl_multifit_fdfsolver_iterate(s);
        if (status)
        {   break; }
        gsl_multifit_gradient(s->J, s->f, g);
        status = gsl_multifit_test_gradient(g, epsabs);
    } while (status == GSL_CONTINUE);
    
    double R = gsl_vector_get(s->x, 0);
    printf("%.12g\n", R);
}


