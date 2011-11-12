#include <config.h>

#include <vector>
#include <string>
#include <sstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#include "coordinate.h"
#include "potential_interpolation.h"

using namespace Cartesian_1d;
namespace c2d = Cartesian_2d;
using namespace Cartesian_3d;

static constexpr char* fn_background_potential = "data/background.potential.out";

Potential_3d* Potential_background(
    const Coordinate_2d& coord, const std::vector<double>& vphi)
{
    size_t n = coord.size();
    assert(n == vphi.size());
    

    auto p2 = std::shared_ptr<Cartesian_2d::Potential_2d>(
        new Cartesian_2d::Potential_interpolation_2d(coord, vphi));
    Potential_3d* p3 = new Potential_rotate(p2);
    return p3;
}

