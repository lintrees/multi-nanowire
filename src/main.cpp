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

struct Nanowire
{
    double x, y;
};

//void Readin_coord(std::vector<Nanowire>& nanowires);
//gsl_matrix* Coord_matrix(const std::vector<Nanowire>& nanowires);
//gsl_matrix* Inner_distance(const gsl_matrix* coord);
gsl_vector* Charge_distri(double phi, double R, const Inner_distance& inner_dist);

int main(int argc, char** argv)
{
    std::ifstream ifs(fn_array_coord);
    Coordinate_2d coord(ifs);
    Inner_distance inner_distance(coord);

//    exit(0);
//
//
//    std::vector<Nanowire> nanowires;
//    Readin_coord(nanowires);
//    gsl_matrix* coord_matrix = Coord_matrix(nanowires);
//    gsl_matrix* inner_distance = Inner_distance(coord_matrix);

    double U0, R;
    assert(argc==3);
    U0 = atof(argv[1]);
    R = atof(argv[2]);
    gsl_vector* Qs = Charge_distri(U0, R, inner_distance);

//    std::cout << coord.size() << std::endl;
//    gsl_vector_fprintf(stdout, Qs, "%g");

//    exit(0);

    size_t n = Qs->size;
    Sptr_Electron_supply pS(new Electron_supply);
    double I_1d_sum;
    std::vector<double> I_1ds(n);
    #pragma omp parallel
    #pragma omp for schedule(dynamic) reduction(+:I_1d_sum)
    for (size_t i = 0; i < n; ++i)
    {
        double Q = gsl_vector_get(Qs, i);
        Sptr_Potential pphi;
        pphi = Sptr_Potential(new Potential_metal_ball(Q, R));
        pphi = Sptr_Potential(new Potential_boost(pphi, -R));
        Sptr_Potential_energy pU(new Potential_energy(pphi, (*pphi)(0)));
        Sptr_Transmission_probability pT(new Transmission_probability(pU));
        Current_density I_1d(pT, pS, Ef, Ef);
        double i1d = I_1d();
        I_1ds[i] = i1d;
        I_1d_sum += i1d;
//        std::cout << phi(0) << "   ..    " << U(2e-9) << std::endl;
//        std::cout << Q << "  " << i1d << std::endl;
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
