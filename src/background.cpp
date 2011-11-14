#include <config.h>

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#include "coordinate.h"
#include "potential_interpolation.h"

using namespace Cartesian_1d;
namespace c2d = Cartesian_2d;
using namespace Cartesian_3d;

static constexpr char fn_background_potential[] = "data/background.potential.out";

Potential_3d* Potential_background()
{
    std::ifstream ifs(fn_background_potential);
    assert(ifs);    
    Coordinate_2d coords(ifs);
    size_t n = coords.size();     
    ifs.clear();
    ifs.seekg(std::ios::beg);
    std::string str;
    std::stringstream ss;
    std::vector<double> vphi;
    double x, y, phi;
    while(getline(ifs, str))
    {
        ss.clear();
        ss.str(str);
        ss >> x >> y >> phi;
        vphi.push_back(phi);
    }
    assert(n == vphi.size());
        auto p2 = std::shared_ptr<Cartesian_2d::Potential_2d>(
        new Cartesian_2d::Potential_interpolation_2d(coords, vphi));
    Potential_3d* p3 = new Potential_rotate(p2);
    return p3;
}

