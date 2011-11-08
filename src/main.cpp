#include <config.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "coordinate.h"
#include "electron_supply.h"
#include "transmission_probability.h"
#include "current.h"

using namespace std;

//using std::cin, std::cout, std::endl;

static const char* fn_array_coord = "data/array_coord.dat";

static const double Ef = -5.;

gsl_vector* Charge_distri(double phi, double R, const Inner_distance& inner_dist);

int main(int argc, char** argv)
{
    std::ifstream ifs(fn_array_coord);
    Coordinate_2d coord(ifs);
    Inner_distance inner_distance(coord);


    double phi0, R;
    assert(argc==3);
    phi0 = atof(argv[1]);
    R = atof(argv[2]);
    gsl_vector* Qs = Charge_distri(phi0, R, inner_distance);

    size_t n = Qs->size;
    Sptr_Electron_supply pS(new Electron_supply);
    double I_1d_sum;
    std::vector<double> I_1ds(n);
    #pragma omp parallel
    #pragma omp for schedule(dynamic) reduction(+:I_1d_sum)
    for (size_t i = 0; i < n; ++i)
    {
        double Q = gsl_vector_get(Qs, i);
        Potential_superimpose<>* pphi = new Potential_superimpose<>;
        Sptr_Potential spphi;        
        spphi = Sptr_Potential(new Potential_metal_sphere(Q, R));
        spphi = Sptr_Potential(new Potential_boost(spphi, -R));
        pphi->push_back(spphi);
        spphi = Sptr_Potential(new Potential_constant(phi0));
        pphi->push_back(spphi);
        double iphi0 = (*pphi)(0);
        spphi = Sptr_Potential(new Potential_metal_sphere_image(R));
        pphi->push_back(spphi);
        spphi = Sptr_Potential(pphi);        
        Sptr_Potential_energy pU(new Potential_energy(spphi, iphi0));
        Sptr_Transmission_probability pT(new Transmission_probability(pU));
        Current_density I_1d(pT, pS, Ef, Ef);
        double i1d = I_1d();
        I_1ds[i] = i1d;
        I_1d_sum += i1d;
    }
    
    std::cout.precision(12);
    for (size_t i = 0; i < n; ++i)
    {
        cout.width(18);
        cout << coord.x(i) << "   ";
        cout.width(18);
        cout << coord.y(i) << "   ";
        cout.width(18);
        cout << gsl_vector_get(Qs, i) << "   ";
        cout.width(18);
        cout  << I_1ds[i]*M_PI*R*R << endl;
    }
    std::cout << I_1d_sum*M_PI*R*R << std::endl;
}
