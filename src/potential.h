#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <cassert>
#include <memory>
#include <vector>

#include <gsl/gsl_math.h>

#define hypot(x,y,z) gsl_hypot3(x,y,z)

class Potential
{
    public:
        virtual double operator() (double x) const  = 0;
};

typedef std::shared_ptr<Potential> Sptr_Potential;

/* metal ball potential */
class Potential_metal_ball: public Potential
{
    public:
        Potential_metal_ball(double Q, double R)
            : _Q(Q), _R(R) { assert(R > 0); }
        double operator() (double x) const;
    private:
        double _Q;
        double _R;
};

/* metal ball potential */
class Potential_metal_ball: public Potential
{
    public:
        Potential_metal_ball(double Q, double R)
            : _Q(Q), _R(R) { assert(R > 0); }
        double operator() (double x) const;
    private:
        double _Q;
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


/*
 * potential in 3d
 */

class Potential_3d
{
    public:
        virtual double operator() (double x, double y, double z) const  {}//= 0;
};
typedef std::shared_ptr<Potential_3d> Sptr_Potential_3d;

/* metal ball potential in 3d*/
class Potential_metal_ball_3d: public Potential_3d
{
    public:
        Potential_metal_ball_3d(double Q, double R)
            : _phi(Q, R) {}   
        double operator() (double x, double y, double z) const
        {
            return _phi(hypot(x, y, z));
        }
    private:    
        Potential_metal_ball _phi;
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

#endif // _POTENTIAL_H_
