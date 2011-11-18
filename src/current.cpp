#include <config.h>

#include <cmath>
#include <cassert>

#include <gsl/gsl_const.h>
#include <gsl/gsl_integration.h>

#include "current.h"

using namespace Cartesian_1d;

static constexpr double epsabs = 0;
static constexpr double epsrel = 1e-4;
static constexpr size_t iter_limit = 200;
static constexpr size_t integration_workspace_size = 4096;
static constexpr int integ_key = GSL_INTEG_GAUSS51;

static constexpr double e = GSL_CONST_MKSA_ELECTRON_CHARGE;


struct int_params_struct{
    const Transmission_probability* pT;
    const Electron_supply* pS;    
    double Ef;
};
static inline double int_func(double E, void* int_params)
{
    int_params_struct* params(static_cast<int_params_struct*>(int_params));
    const Transmission_probability* pT = params->pT;    
    const Electron_supply* pS = params->pS;
    const double Ef = params->Ef;
    return (*pS)(E-Ef) * (*pT)(E);
}
double Current_density::operator() () const
{
    gsl_integration_workspace* w
        = gsl_integration_workspace_alloc(integration_workspace_size);
    double int_Emax = _Ef + _pS->Emax();
    int_params_struct int_params = {_pT.get(), _pS.get(), _Ef};
    gsl_function f;
    f.function = &int_func;
    f.params = &int_params;    
    double result, abserr;
    gsl_integration_qag(&f, _int_Emin, int_Emax, epsabs, epsrel, iter_limit,
        integ_key, w, &result, &abserr);
    gsl_integration_workspace_free(w);
    return e*e*result;
}

