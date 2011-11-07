#include <config.h>

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>

#include <gsl/gsl_const.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "potential.h"

static const size_t N = 200;
static const double xmin = -1.5e-6, xmax = 1.5e-6;
static const double zmin = 0, zmax = 1e-6;


struct Nanowire
{
    double x, y;  
};

void Readin_coord(std::vector<Nanowire>& nanowires);
gsl_matrix* Coord_matrix(const std::vector<Nanowire>& nanowires);
gsl_matrix* Inner_distance(const gsl_matrix* coord);
gsl_vector* Charge_distri(double phi, double R, const gsl_matrix* inner_dist);

int main(int argc, char** argv)
{
    std::vector<Nanowire> nanowires;
    Readin_coord(nanowires);
    gsl_matrix* coord = Coord_matrix(nanowires);
    gsl_matrix* inner_distance = Inner_distance(coord);    
    double U0, R;
    assert(argc==3);
    U0 = atof(argv[1]);
    R = atof(argv[2]);    
    gsl_vector* Qs = Charge_distri(U0, R, inner_distance);
    
    
    size_t n = Qs->size;
    Potential_superimpose_3d<> phi;
    for (size_t i = 0; i < n; ++i)
    {
        double Q = gsl_vector_get(Qs, i);
        double x = gsl_matrix_get(coord, i, 0);
        double y = gsl_matrix_get(coord, i, 1);
        Sptr_Potential_3d pphi_i;
        pphi_i = Sptr_Potential_3d(new Potential_metal_ball_3d(Q, R));
        pphi_i = Sptr_Potential_3d(new Potential_boost_3d(pphi_i, x, y, 0));
        phi.push_back(pphi_i);
    }
    
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
