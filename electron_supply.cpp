#include <config.h>

#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include <gsl/gsl_const.h>
#include <gsl/gsl_spline.h>


static const char* fn_fztbl = "fztbl.dat";
static const double start = -.5;
static double stop;
static const double step = 1e-3;

static gsl_spline* spline;
static gsl_interp_accel* acc;

void electron_supply_init()
{
    std::ifstream ifs(fn_fztbl);
    assert(ifs);
    double x = start;
    double y;
    std::vector<double> xi, yi;
    while (ifs >> y)    
    {   
        xi.push_back(x);
        yi.push_back(y);
        x += step;        
    }
    stop = x;
    size_t n = xi.size();
    spline = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(spline, xi.data(), yi.data(), n);
    acc = gsl_interp_accel_alloc();
}

double electron_supply(double x)
{
    assert(start <= x && x <= stop);
    return gsl_spline_eval(spline, x, acc);
}

int main()
{
    electron_supply_init();
    std::cout << electron_supply(-0.40056) << std::endl;
}
