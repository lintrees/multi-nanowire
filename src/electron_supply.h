#ifndef _ELECTRON_SUPPLY_H_
#define _ELECTRON_SUPPLY_H_

#include <memory>

#include <gsl/gsl_spline.h>

class Electron_supply
{
    public:
        Electron_supply();
        ~Electron_supply();
        double operator() (double E) const;
        double Emin() const { return _Emin; }
        double Emax() const { return _Emax; }        
    private:
        gsl_spline* _spline;
        gsl_interp_accel* _acc;
        double _Emin, _Emax;
};
typedef std::shared_ptr<Electron_supply> Sptr_Electron_supply;

#endif // _ELECTRON_SUPPLY_H_
