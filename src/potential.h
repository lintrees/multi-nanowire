#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <cassert>
#include <vector>

#include <gsl/gsl_math.h>

#define hypot(x,y,z) gsl_hypot3(x,y,z)


class Potential
{
    public:
        virtual double operator() (double x)  const {}// = 0;
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
        typedef Potential_boost Base;
        Potential_boost(const Potential& phi, double x0)
            : _phi(phi), _x0(x0) {}
        double operator() (double x) const
        {   return _phi(x-_x0); }
    private:
        const Potential _phi;
        double _x0;
};


/*
 * potential in 3d
 */

class Potential_3d
{
    public:
        virtual double operator() (double x, double y, double z) const  {}//= 0;
};

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
        Potential_boost_3d(const Potential_3d& phi, double x0, double y0, double z0)
            : _phi(phi), _x0(x0), _y0(y0), _z0(z0) {}
        double operator()(double x, double y, double z) const
        {
            return _phi(x-_x0, y-_y0, z-_z0);
        }
    private:
        Potential_3d _phi;
        double _x0, _y0, _z0;    
};

/*  superimpose potential in 3d */
/*
template <class Container = std::vector<const Potential_3d>>
class Potential_superimpose_3d: public Potential_3d, public Container
{
    public:
        double operator()(double x, double y, double z) const
        {
            double sum = 0;
            for (auto iPhi = this->cbegin(); iPhi != this->cend(); ++iPhi)
            {   sum += (*iPhi)(x, y, z); }
            return sum;
        }
};
*/
class Potential_energy
{
    public:
        Potential_energy(const Potential& phi, double phi0)
            : _phi(phi), _phi0(phi0) {}
        double operator() (double x) const
        {
            return - (_phi(x)-_phi0);
        }
    private:
        Potential _phi;
        double _phi0;
};
#endif // _POTENTIAL_H_
