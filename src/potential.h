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
        virtual double operator() (double x) const;
    private:
        
};

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

class Potential_boost: public Potential
{
    public:
        Potential_boost(const Potential& Phi, double x0): _spPhi(new Phi), _x0(x0) {}
        double operator() (double x) const
        {
            return (*_spPhi)(x-_x0);
        }
    private:
        Shared_ptr_Potential _spPhi;
        double _x0;
};


class Potential_3d
{
    public:
        virtual double operator() (double x, double y, double z) const = 0;
};





class Potential_metal_ball_3d: public Potential_3d
{
    public:
        Potential_metal_ball_3d(double Q, double R) : _phi(Q, R) {}
        double operator() (double x, double y, double z) const
        {
            return _phi(hypot(x, y, z));
        }
    private:
        Potential_metal_ball _phi;
};


class Potential_boost_3d: public Potential_3d
{
    public:
        Potential_boost_3d(const Potential_3d& Phi, double x0, double y0, double z0)
            : _pPhi(&Phi), _x0(x0), _y0(y0), _z0(z0) {}
        double operator()(double x, double y, double z) const
        {
            return (*_pPhi)(x-_x0, y-_y0, z-_z0);
        }
    private:
        const Potential_3d* _pPhi;
        double _x0, _y0, _z0;
};

template <class Container = std::vector<const Potential_3d*> >
class Potential_STL_3d: public Potential_3d, public Container
{
    public:
        double operator()(double x, double y, double z) const
        {
            double sum = 0;
            for (auto ipPhi = this->begin(); ipPhi != this->end(); ++ipPhi)
            {   sum += (**ipPhi)(x, y, z); }
            return sum;
        }
};

class Potential_energy
{
    public:
        Potential_energy(const Potential& Phi, double Ec): _pPhi(&Phi), _Ec(Ec) {}
        double operator() (double x) const
        {
            return _Ec - (*_pPhi)(x);
        }
    private:
        const Potential* _pPhi;
        double _Ec;
};

#endif // _POTENTIAL_H_
