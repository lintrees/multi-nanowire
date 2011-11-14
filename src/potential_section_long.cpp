#include <config.h>

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>

#include <gsl/gsl_const.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "coordinate.h"
#include "potential.h"

static const char* fn_array_coord = "data/array_coord.dat";

static const size_t N = 200;
static const double xmin = -1.3e-6, xmax = 1.3e-6;
static const double zmin = 0, zmax = .9e-6;

using namespace std;
using namespace Cartesian_1d;
using namespace Cartesian_3d;

static constexpr double hNW = 1e-6;

typedef std::shared_ptr<Potential_3d> Sp3d;
typedef std::shared_ptr<Potential> Sp1d;
typedef std::shared_ptr<Potential_superimpose<>> Sp1d_sup;
typedef std::shared_ptr<Potential_superimpose_3d<>> Sp3d_sup;

vector<double> Charge_distri(const vector<double>& vphi, double R, const Inner_distance& inner_dist);
Potential_3d* Potential_background();

int main(int argc, char** argv)
{
    double phi0, R;
    assert(argc==3);
    phi0 = atof(argv[1]);
    R = atof(argv[2]);
    std::ifstream ifs(fn_array_coord);
    Coordinate_2d coord(ifs);
    ifs.close();
    size_t n = coord.size();
    Inner_distance inner_distance(coord);
    Potential_3d* pphi_background = Potential_background();
    Sp3d spphi_background(Potential_background());
    vector<double> vphi_background(n), vphi_delta(n);
    for (int i = 0; i < n; ++i)
    {
        vphi_background[i] = (*spphi_background)(coord.x(i), coord.y(i), hNW);
        vphi_delta[i] = phi0 - vphi_background[i];
    }
    vector<double> vQ = Charge_distri(vphi_delta, R, inner_distance);
    
//    std::shared_ptr<Potential_superimpose_3d<>>
//        pphi_sup(new Potential_superimpose_3d<>);
//    for (size_t i = 0; i < n; ++i)
//    {
//        double Q = vQ[i];
//        double x = coord.x(i);
//        double y = coord.y(i);
//        Sptr_Potential_3d pphi_i;
//        pphi_i = Sptr_Potential_3d(new Potential_metal_sphere_3d(Q, R));
//        pphi_i = Sptr_Potential_3d(new Potential_boost_3d(pphi_i, x, y, 0));
//        pphi_sup->push_back(pphi_i);
//    }
//    Potential_shift_3d phi(pphi_sup, phi0);

    Sp3d_sup spphi_sup(new Potential_superimpose_3d<>);
    // add all metal sphere potential
    for (size_t i = 0; i < n; ++i)
    {
        double Q = vQ[i];
        double x = coord.x(i);
        double y = coord.y(i);
        Sp3d spphi;
        spphi = Sp3d(new Potential_metal_sphere_3d(Q, R));
        spphi = Sp3d(new Potential_boost_3d(spphi, x, y, 0));
        spphi_sup->push_back(spphi);
    }
    // add all background potential
    Sp3d spphi_bg_boost(new Potential_boost_3d(spphi_background, 0, 0, -hNW));
    spphi_sup->push_back(spphi_bg_boost);    
    
    double dx = (xmax-xmin)/N;
    double dz = (zmax-zmin)/N;
    for (double x = xmin; x <= xmax; x += dx)
    {   for (double z = zmin; z <= zmax; z += dz)
        {
            double p = (*spphi_sup)(x, 0, z);
            std::cout << x << "   "
                      << z << "   "
                      << p << std::endl;
        }
    }

    return 0;
}
