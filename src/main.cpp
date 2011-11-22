#include <config.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "coordinate.h"
#include "electron_supply.h"
#include "transmission_probability.h"
#include "current.h"

using namespace std;
using namespace Cartesian_1d;
using namespace Cartesian_3d;

static constexpr char fn_array_coord[] = "data/array_coord.dat";
static constexpr char fn_background_potential[] = "data/background.potential.out";
static constexpr double Ef = -5.;
static constexpr double phi_0 =(-5)-(-4.05); // zero point

typedef std::shared_ptr<Potential_3d> Sp3d;
typedef std::shared_ptr<Potential> Sp1d;
typedef std::shared_ptr<Potential_superimpose<>> Sp1d_sup;
typedef std::shared_ptr<Potential_superimpose_3d<>> Sp3d_sup;

vector<double> Charge_distri(const vector<double>& vphi, double R, const Inner_distance& inner_dist);
Potential_3d* Potential_background(std::istream& is);

int main(int argc, char** argv)
{
    double phi0, R, hNW;
    double Ec; // the bottom of the conduction band
    assert(argc==5);
    phi0 = atof(argv[1]);
    R = atof(argv[2]);
    Ec = atof(argv[3]);
    hNW = atof(argv[4]);

    std::ifstream ifs(fn_array_coord);
    Coordinate_2d coord(ifs);
    ifs.close();
    size_t n = coord.size();
    Inner_distance inner_distance(coord);
    
    ifs.open(fn_background_potential);
    Sp3d spphi_background(Potential_background(ifs));
    ifs.close();
    vector<double> vphi_background(n), vphi_delta(n);
    for (int i = 0; i < n; ++i)
    {
        vphi_background[i] = (*spphi_background)(coord.x(i), coord.y(i), hNW);
        vphi_delta[i] = phi0 - vphi_background[i];
    }
    vector<double> vQ = Charge_distri(vphi_delta, R, inner_distance);

    Sptr_Electron_supply pS(new Electron_supply);
    double I_1d_sum(0);

#ifndef SINGLE_SPHERE_CHARGE_APPROX
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
#endif // SINGLE_SPHERE_CHARGE_APPROX

    vector<double> vI_1d(n);
    vector<shared_ptr<Potential_energy> > vU(n);
    #pragma omp parallel
    #pragma omp for schedule(dynamic) reduction(+:I_1d_sum)
    for (size_t i = 0; i < n; ++i)
    {
        double Q = vQ[i];
        double x = coord.x(i);
        double y = coord.y(i);
        Sp1d_sup spphi_sup_i(new Potential_superimpose<>);

#ifdef SINGLE_SPHERE_CHARGE_APPROX
        Sp1d spphi;
        spphi = Sp1d(new Potential_metal_sphere(Q, R));
        spphi = Sp1d(new Potential_boost(spphi, -R));
        spphi_sup_i->push_back(spphi);
#else
        spphi_sup_i->push_back(Sp1d(new Potential_path(spphi_sup, x, y, R, 0, 0)));
#endif // SINGLE_SPHERE_CHARGE_APPROX

        double phi0i = (*spphi_sup_i)(0);
        spphi_sup_i->push_back(Sp1d(new Potential_metal_sphere_image(R)));
        shared_ptr<Potential_energy> pU(new Potential_energy(spphi_sup_i, phi0i+phi_0));
        vU[i] = pU;
        shared_ptr<Transmission_probability> pT(new Transmission_probability(pU));
        Current_density I_1d(pT, pS, Ef, Ec);
        double i1d = I_1d();
        vI_1d[i] = i1d;
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
        cout << vQ[i] << "   ";
        cout.width(18);
        cout << vI_1d[i]*M_PI*R*R << endl;
    }
    std::cout << I_1d_sum*M_PI*R*R << std::endl;
    
#ifndef NDEBUG
    size_t imax;
    double maxI_1d = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (vI_1d[i] > maxI_1d)
        {   imax = i; }
    }    
    for (double x = 0; x <= 5e-9; x+=0.1e-9)
    {
        std::cerr << x << "    " << (*vU[imax])(x) << endl;
    }
#endif // NDEBUG
}
