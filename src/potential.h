#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <cassert>
#include <memory>
#include <vector>

#include <gsl/gsl_math.h>

//#define hypot(x,y,z) gsl_hypot3(x,y,z)

/*
 * Potential in 1d Cartesian
 */
namespace Cartesian_1d
{
class Potential
{
    public:
        virtual double operator() (double x) const  = 0;
};

typedef std::shared_ptr<Potential> Sptr_Potential;

/* constant potential */
class Potential_constant: public Potential
{
    public:
        Potential_constant(double phi0)
            : _phi0(phi0) {}
        double operator() (double x) const
        {   return _phi0; }
    private:
        double _phi0;
};

/* metal sphere potential */
class Potential_metal_sphere: public Potential
{
    public:
        Potential_metal_sphere(double Q, double R)
            : _Q(Q), _R(R) { assert(R > 0); }
        double operator() (double x) const;
        /* distance from the centre */
    private:
        double _Q;
        double _R;
};

/* metal sphere image potential */
class Potential_metal_sphere_image : public Potential
{
    public:
        Potential_metal_sphere_image(double R)
            : _R(R) { assert(R > 0); }
        double operator() (double d) const;
        /* distance from the sphere */
    private:
        double _R;
};

/* boost potential */
class Potential_boost: public Potential
{
    public:
        Potential_boost(const Sptr_Potential& pphi, double x0)
            : _pphi(pphi), _x0(x0) {}
        double operator() (double x) const
        {   return (*_pphi)(x-_x0); }
    private:
        const Sptr_Potential _pphi;
        double _x0;
};

/* shift potential */
class Potential_shift: public Potential
{
    public:
        Potential_shift(const Sptr_Potential& pphi, double dphi)
            : _pphi(pphi), _dphi(dphi) {}
        double operator() (double x) const
        {   return (*_pphi)(x) + _dphi; }
    private:
        const Sptr_Potential _pphi;
        double _dphi;
};

/*  superimpose potential */
template <class Container = std::vector<Sptr_Potential>>
class Potential_superimpose: public Potential, public Container
{
    public:
        double operator()(double x) const
        {
            double sum = 0;
            for (auto ipPhi = this->cbegin(); ipPhi != this->cend(); ++ipPhi)
            {   sum += (**ipPhi)(x); }
            return sum;
        }
};

} //namespace Cartesian_1d
/*
 * potential in 3d
 */
namespace Cartesian_3d
{

class Potential_3d
{
    public:
        virtual double operator() (double x, double y, double z) const = 0;
};
typedef std::shared_ptr<Potential_3d> Sptr_Potential_3d;

/* metal sphere potential in 3d*/
class Potential_metal_sphere_3d: public Potential_3d
{
    public:
        Potential_metal_sphere_3d(double Q, double R)
            : _phi(Q, R) {}
        double operator() (double x, double y, double z) const
        {
            return _phi(gsl_hypot3(x, y, z));
        }
    private:
        const Cartesian_1d::Potential_metal_sphere _phi;
};

/*  boost potential in 3d */
class Potential_boost_3d: public Potential_3d
{
    public:
        Potential_boost_3d(const Sptr_Potential_3d& pphi, double x0, double y0, double z0)
            : _pphi(pphi), _x0(x0), _y0(y0), _z0(z0) {}
        double operator()(double x, double y, double z) const
        {
            return (*_pphi)(x-_x0, y-_y0, z-_z0);
        }
    private:
        const Sptr_Potential_3d _pphi;
        double _x0, _y0, _z0;
};

/* shift potential in 3d */
class Potential_shift_3d: public Potential_3d
{
    public:
        Potential_shift_3d(const Sptr_Potential_3d& pphi, double dphi)
            : _pphi(pphi), _dphi(dphi) {}
        double operator()(double x, double y, double z) const
        {   return (*_pphi)(x, y, z) + _dphi; }
    private:
        const Sptr_Potential_3d _pphi;
        double _dphi;
};

/* reflect to xy-plane potential in 3d */
class Potential_reflect_xy_3d: public Potential_3d
{
    public:
        Potential_reflect_xy_3d(const Sptr_Potential_3d& pphi, double z0=0)
            : _pphi(pphi), _z0_2(z0*2) {}
        double operator()(double x, double y, double z) const
        {   return (*_pphi)(x, y, _z0_2-z); }
    private:
        const Sptr_Potential_3d _pphi;
        double _z0_2;
};

/* opposite potential in 3d */
class Potential_opposite_3d: public Potential_3d
{
    public:
        Potential_opposite_3d(const Sptr_Potential_3d& pphi)
            : _pphi(pphi) {}
        double operator()(double x, double y, double z) const
        {   return - (*_pphi)(x, y, z); }
    private:
        const Sptr_Potential_3d _pphi;
};

/*  superimpose potential in 3d */
template <class Container = std::vector<Sptr_Potential_3d>>
class Potential_superimpose_3d: public Potential_3d, public Container
{
    public:
        double operator()(double x, double y, double z) const
        {
            double sum = 0;
            for (auto ipPhi = this->cbegin(); ipPhi != this->cend(); ++ipPhi)
            {   sum += (**ipPhi)(x, y, z); }
            return sum;
        }
};

} // namespace Cartesian_3d


namespace Cartesian_1d
{
/* 1d path potential from 3d */
class Potential_path: public Potential
{
    public:
        Potential_path(const Cartesian_3d::Sptr_Potential_3d& pphi,
            double x0, double y0, double z0, double theta, double phi)
            : _pphi(pphi), _x0(x0), _y0(y0), _z0(z0),
                _dx(sin(theta)*cos(phi)), _dy(sin(theta)*sin(phi)), _dz(cos(theta)) {}
        double operator() (double r) const
        {   return (*_pphi)(r*_dx+_x0, r*_dy+_y0, r*_dz+_z0); }
    private:
        const Cartesian_3d::Sptr_Potential_3d _pphi;
        double _x0, _y0, _z0;
        double _dx, _dy, _dz;
};


/* Potential energy */
class Potential_energy
{
    public:
        Potential_energy(const Sptr_Potential& pphi, double phi0)
            : _pphi(pphi), _phi0(phi0) {}
        double operator() (double x) const
        {
            return - ((*_pphi)(x)-_phi0);
        }
    private:
        const Sptr_Potential _pphi;
        double _phi0;
};
typedef std::shared_ptr<Potential_energy> Sptr_Potential_energy;

} // namespace Cartesian_1d

#endif // _POTENTIAL_H_
