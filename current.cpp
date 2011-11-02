#include <config.h>

#include <cmath>
#include <cassert>

#include <gsl/gsl_const.h>


static const double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
static const double epsilon_0_4_pi = epsilon_0 * 4*M_PI;


static inline double potential(double x, double Q, double R)
/* electric charge, ball radius */
{
    assert(x >= R);
    return Q / (epsilon_0_4_pi * (x - R));
}


