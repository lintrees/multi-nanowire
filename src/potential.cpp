#include <config.h>

#include <cmath>
#include <cassert>

#include <gsl/gsl_const.h>

#include "potential.h"

static const double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
static const double epsilon_0_4_pi = epsilon_0 * 4*M_PI;


double Potential_base_metal_ball::operator() (double x) const
/* distance from origin */
{
    assert(x >= 0);
    if (x > _R)
    {   return _Q / (epsilon_0_4_pi * x); }
    else
    {   return _Q / (epsilon_0_4_pi * _R); }
}

