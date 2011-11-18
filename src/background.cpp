#include <config.h>

#include <vector>
#include <string>
#include <sstream>
#include <istream>

#include "coordinate.h"
#include "potential_interpolation.h"

namespace c1d = Cartesian_1d;
namespace c2d = Cartesian_2d;
namespace c3d = Cartesian_3d;

c3d::Potential_3d* Potential_background(std::istream& is)
{
    Coordinate_2d coords(is);
    size_t n = coords.size();  
    is.clear();   
    is.seekg(std::ios::beg);
    std::string str;
    std::stringstream ss;
    std::vector<double> vphi;
    double x, y, phi;
    while(getline(is, str))
    {
        ss.clear();
        ss.str(str);
        ss >> x >> y >> phi;
        vphi.push_back(phi);
    }
    assert(vphi.size() == n);
    std::shared_ptr<c2d::Potential_2d>
        p2(new c2d::Potential_interpolation_2d(coords, vphi));
    c3d::Potential_3d* p3 = new c3d::Potential_rotate(p2);
    return p3;
}

