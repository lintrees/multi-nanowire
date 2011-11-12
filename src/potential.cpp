#include <config.h>

#include <cmath>
#include <cassert>

#include <gsl/gsl_const.h>

#include "potential.h"

static constexpr double eps_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
static constexpr double e = GSL_CONST_MKSA_ELECTRON_CHARGE;


namespace Cartesian_1d
{
double Potential_metal_sphere::operator() (double x) const
/* distance from the centre */
{
    constexpr static double coeff = 1/(eps_0 *4*M_PI);
    assert(x >= 0);
    if (x > _R)
    {   return coeff* _Q/x; }
    else
    {   return coeff* _Q/_R; }
}

double Potential_metal_sphere_image::operator() (double d) const
/* distance from the sphere */
{
    constexpr static double coeff = +e/(eps_0 *4*M_PI);
    return coeff* _R/(2*d*(d+2*_R));
}
} // namespace Cartesian_1d
