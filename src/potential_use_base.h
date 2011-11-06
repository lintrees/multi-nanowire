#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <cassert>
#include <vector>
#include <memory>

#include <gsl/gsl_math.h>

#define hypot(x,y,z) gsl_hypot3(x,y,z)


class Potential_base
{
    public:
        virtual double operator() (double x) const = 0;
};

typedef std::shared_ptr<Potential_base> Shared_ptr_Potential_base;

class Potential
{
    public:
        virtual double operator() (double x) const {return -1e100;}
//        {
//            return (*_pbase)(x);
//        }
//    public:
//        Shared_ptr_Potential_base _pbase;
};

/* metal ball potential */
class Potential_base_metal_ball: public Potential_base
{
    public:
        Potential_base_metal_ball(double Q, double R)
            : _Q(Q), _R(R) { assert(R > 0); }
        double operator() (double x) const;
    private:
        double _Q;
        double _R;
};
class Potential_metal_ball: public Potential
{
    public:
        typedef Potential_base_metal_ball Base;
        Potential_metal_ball(double Q, double R)
            : _pbase(new Base(Q, R)) {}
        double operator() (double x) const
        {   return (*_pbase)(x); }
    private:
        std::shared_ptr<Base> _pbase;
};

/* boost potential */
class Potential_base_boost: public Potential_base
{
    public:
        Potential_base_boost(const Potential& phi, double x0)
            : _pphi(new Potential(phi)), _x0(x0) {}
        double operator() (double x) const
        {
            return (*_pphi)(x-_x0);
        }
    private:
        Potential* _pphi;
        double _x0;
};
class Potential_boost: public Potential
{
    public:
        typedef Potential_base_boost Base;
        Potential_boost(const Potential& phi, double x0)
            : _pbase(new Base(phi, x0)) {}
        double operator() (double x) const
        {   return (*_pbase)(x); }
    private:
        std::shared_ptr<Base> _pbase;
};


/*
 * potential in 3d
 */

class Potential_base_3d
{
    public:
        virtual double operator() (double x, double y, double z) const = 0;
};

typedef std::shared_ptr<Potential_base_3d> Shared_ptr_Potential_base_3d;

class Potential_3d
{
    public:
        double operator() (double x, double y, double z) const
        {
            return (*_pbase)(x, y, z);
        }
    protected:
        Shared_ptr_Potential_base_3d _pbase;
};

/* metal ball potential in 3d*/
class Potential_base_metal_ball_3d: public Potential_base_3d
{
    public:
        Potential_base_metal_ball_3d(double Q, double R)
            : _phi(Q, R) {}
        double operator() (double x, double y, double z) const
        {
            return _phi(hypot(x, y, z));
        }
    private:
        Potential_metal_ball _phi;
};
class Potential_metal_ball_3d: public Potential_3d
{
    public:
        Potential_metal_ball_3d(double Q, double R)
            : _pbase(new Potential_base_metal_ball_3d(Q, R)) {}
        
    private:
        Shared_ptr_Potential_base_3d _pbase;
};

/*  boost potential in 3d */
class Potential_base_boost_3d: public Potential_base_3d
{
    public:
        Potential_base_boost_3d(const Potential_3d& phi, double x0, double y0, double z0)
            : _phi(phi), _x0(x0), _y0(y0), _z0(z0) {}
        double operator()(double x, double y, double z) const
        {
            return _phi(x-_x0, y-_y0, z-_z0);
        }
    private:
        const Potential_3d _phi;
        double _x0, _y0, _z0;
};
class Potential_boost_3d: public Potential_3d
{
    public:
        Potential_boost_3d(const Potential_3d& phi, double x0, double y0, double z0)
            : _pbase(new Potential_base_boost_3d(phi, x0, y0, z0)) {}
    private:
        Shared_ptr_Potential_base_3d _pbase;        
};

/*  superimpose potential in 3d */
template <class Container>
class Potential_base_superimpose_3d: public Potential_base_3d, public Container
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
template <class Container = std::vector<const Potential_3d> >
class Potential_superimpose_3d: public Potential_3d
{
    public:
        typedef Potential_base_superimpose_3d<Container> Base;
//        Potential_superimpose_3d(): _pbase(new Base) {}
//        void push_back(const Potential_3d& phi)
//        {
//            _pbase->push_back(phi);
//        }
    private:
        std::shared_ptr<Base> _pbase;
};

class Potential_energy_base
{
    public:
        Potential_energy_base(const Potential& phi, double phi0): _pphi(new Potential(phi)), _phi0(phi0) {}
        double operator() (double x) const
        {
            return - ((*_pphi)(x)-_phi0);
        }
    private:
        const Potential* _pphi;
        double _phi0;
};
class Potential_energy
{
    public:
        typedef Potential_energy_base Base;
        Potential_energy(const Potential& phi, double phi0)
            : _pbase(new Base(phi, phi0)) {}
        double operator() (double x) const
        {
            return (*_pbase)(x);
        }
    private:
        std::shared_ptr<Base> _pbase;
};
#endif // _POTENTIAL_H_
