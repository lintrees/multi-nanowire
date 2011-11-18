#include <config.h>

#include <iostream>
#include <fstream>

#include "potential.h"


namespace c3d = Cartesian_3d;
namespace c1d = Cartesian_1d;

typedef std::shared_ptr<c3d::Potential_3d> Sp3d;
typedef std::shared_ptr<c1d::Potential> Sp1d;

static constexpr char fn_single_potential[] = "data/single.potential.out";
static constexpr double length = 30e-9;
static constexpr size_t n = 100;

c3d::Potential_3d* Potential_background(std::istream& is);

int main(int argc, char** argv)
{
    assert(argc==3);
    double x0, y0;
    x0 = atof(argv[1]);
    y0 = atof(argv[2]);

    std::ifstream ifs(fn_single_potential);
    Sp3d pphi_bg(Potential_background(ifs));
    ifs.close();

    Sp1d pphi_path(new c1d::Potential_path(pphi_bg, x0, 0, y0, 0, 0));

    double dx = length/n;
    for (double x = 0; x < length; x+= dx)
    {
        std::cout << x << "    " << (*pphi_path)(x) << std::endl;
    }
    return 0;
}
