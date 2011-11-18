#include <config.h>

#include <cmath>

#include <gsl/gsl_const.h>
#include <gsl/gsl_integration.h>

#include "transmission_probability.h"

using namespace Cartesian_1d;

static constexpr double xmax = 5e-9; // max potential barriar length
static constexpr double epsabs = 0;
static constexpr double epsrel = 1e-8;
static constexpr size_t iter_limit = 100;
static constexpr size_t integration_workspace_size = 1024;
//static constexpr int integ_key = GSL_INTEG_GAUSS41;

static constexpr double m_e = GSL_CONST_MKSA_MASS_ELECTRON;
static constexpr double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
static constexpr double eV = GSL_CONST_MKSA_ELECTRON_VOLT;


struct int_params_struct{
    const Potential_energy* pU;
    double E;
};
static inline double int_func(double x, void* int_params)
{
    static const double coeff = eV*2*m_e/(hbar*hbar);
    const Potential_energy* pU;
    pU = static_cast<int_params_struct*>(int_params)->pU;
    double E = static_cast<int_params_struct*>(int_params)->E;
    double U_x = (*pU)(x);
    if (E >= U_x)
    {   return 0; }
    else
    {   return sqrt(coeff*(U_x-E)); }
}
double Transmission_probability::operator() (double E) const
/*  T(E) = \exp{-2\int \d x \sqrt{2*m_e/\hbar^2 *(U(x)-E)}}
    E in (eV)
*/
{
    gsl_integration_workspace* w;
    w = gsl_integration_workspace_alloc(integration_workspace_size);
    int_params_struct int_params = {_pU.get(), E};
    gsl_function f;
    f.function = &int_func;
    f.params = &int_params;
    double result, abserr;
    gsl_integration_qags(&f, 0, xmax, epsabs, epsrel, iter_limit,
        /*integ_key,*/ w, &result, &abserr);
    gsl_integration_workspace_free(w);
    if (result == 0)
    {  return 1; }
    else
    {  return exp(-2*result);}
}
