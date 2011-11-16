#include <config.h>

#include <omp.h>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>

#include "coordinate.h"
#include "potential.h"

static constexpr char fn_array_coord[] = "data/array_coord.dat";

static constexpr size_t N = 200;
static constexpr double xmin = -1.4e-6, xmax = 1.4e-6;
static constexpr double zmin = 0, zmax = 2e-6;

static constexpr double hNW = 1e-6;

using namespace std;
using namespace Cartesian_1d;
using namespace Cartesian_3d;

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
    Sp3d spphi_background(Potential_background());
    vector<double> vphi_background(n), vphi_delta(n);
    for (int i = 0; i < n; ++i)
    {
        vphi_background[i] = (*spphi_background)(coord.x(i), coord.y(i), hNW);
        vphi_delta[i] = phi0 - vphi_background[i];
    }
    vector<double> vQ = Charge_distri(vphi_delta, R, inner_distance);

    Sp3d_sup spphi_sup(new Potential_superimpose_3d<>);
    // add all metal sphere potential
    Sp3d_sup spphi_sup_sphere(new Potential_superimpose_3d<>);
    for (size_t i = 0; i < n; ++i)
    {
        double Q = vQ[i];
        double x = coord.x(i);
        double y = coord.y(i);
        Sp3d spphi(new Potential_metal_sphere_3d(Q, R));
        spphi = Sp3d(new Potential_boost_3d(spphi, x, y, hNW));
        spphi_sup_sphere->push_back(spphi);
    }
    spphi_sup->push_back(spphi_sup_sphere);
    // add induced charge potential
    Sp3d spphi_refl(new Potential_reflect_xy_3d(spphi_sup_sphere, 0));
    Sp3d spphi_opp(new Potential_opposite_3d(spphi_refl));
    spphi_sup->push_back(spphi_opp);        
    // add all background potential
//    Sp3d spphi_bg_boost(new Potential_boost_3d(spphi_background, 0, 0, -hNW));
    spphi_sup->push_back(spphi_background);    
    
    double dx = (xmax-xmin)/N;
    double dz = (zmax-zmin)/N;
    double grid[N+1][N+1];
    #pragma omp parallel
    for (int i = 0; i <= N; ++i)
    {
        double x = xmin + i*dx;
        #pragma omp for schedule(guided) nowait
        for (int j = 0; j <= N; ++j)
        {
            double z = zmin +j*dz;
            grid[i][j] = (*spphi_sup)(x, 0, z);
        }
    }
    #pragma omp barrier
    
    #pragma omp single
    for (int i = 0; i <= N; ++i)
    {
        double x = xmin + i*dx;
        for (int j = 0; j <= N; ++j)
        {
            double z = zmin +j*dz;
            std::cout << x << "   "
                      << z << "   "
                      << grid[i][j] << std::endl;
        }
    }

    return 0;
}
