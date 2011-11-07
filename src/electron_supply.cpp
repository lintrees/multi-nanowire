#include <config.h>

#include <cassert>
#include <fstream>
#include <vector>

#include <gsl/gsl_spline.h>

#include "electron_supply.h"


static const char* fn_fztbl = "data/fztbl.dat";
static const double start = -.5;
static const double step = 1e-3;

Electron_supply::Electron_supply()
{
    std::ifstream ifs(fn_fztbl);
    assert(ifs);
    double x = start;
    double y;
    std::vector<double> X, Y;
    while (ifs >> y)    
    {   
        X.push_back(x);
        Y.push_back(y);
        x += step;
    }
    _Emin = start;
    _Emax = x - step;
    size_t n = X.size();
    _spline = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(_spline, X.data(), Y.data(), n);
    _acc = gsl_interp_accel_alloc();
}
Electron_supply::~Electron_supply()
{
    gsl_interp_accel_free(_acc);
    gsl_spline_free(_spline);
}
double Electron_supply::operator() (double E) const
{
    assert(_Emin <= E && E <= _Emax);
    return gsl_spline_eval(_spline, E, _acc);
}
