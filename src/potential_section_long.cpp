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
static const double xmin = -1.5e-6, xmax = 1.5e-6;
static const double zmin = 0, zmax = 1e-6;


struct Nanowire
{
    double x, y;  
};

gsl_vector* Charge_distri(double phi, double R, const Inner_distance& inner_dist);

int main(int argc, char** argv)
{

    std::ifstream ifs(fn_array_coord);
    Coordinate_2d coord(ifs);
    Inner_distance inner_distance(coord);
    double U0, R;
    assert(argc==3);
    U0 = atof(argv[1]);
    R = atof(argv[2]);
    gsl_vector* Qs = Charge_distri(U0, R, inner_distance);
    
    size_t n = Qs->size;
    std::shared_ptr<Potential_superimpose_3d<>>
        pphi_sup(new Potential_superimpose_3d<>);
    for (size_t i = 0; i < n; ++i)
    {
        double Q = gsl_vector_get(Qs, i);
        double x = coord.x(i);
        double y = coord.y(i);
        Sptr_Potential_3d pphi_i;
        pphi_i = Sptr_Potential_3d(new Potential_metal_sphere_3d(Q, R));
        pphi_i = Sptr_Potential_3d(new Potential_boost_3d(pphi_i, x, y, 0));
        pphi_sup->push_back(pphi_i);
    }
    Potential_shift_3d phi(pphi_sup, U0);        
    
    double dx = (xmax-xmin)/N;
    double dz = (zmax-zmin)/N;

    for (double x = xmin; x <= xmax; x += dx)
    {   for (double z = zmin; z <= zmax; z += dz)
        {
            double p = phi(x, 0, z);
            std::cout << x << "   "
                      << z << "   "
                      << p << std::endl;
        }
    }

    return 0;
}
