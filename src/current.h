#ifndef _CURRENT_H_
#define _CURRENT_H_

#include "transmission_probability.h"
#include "electron_supply.h"


class Current_density
{
    public:
        Current_density(const Sptr_Transmission_probability pT,
            const Sptr_Electron_supply pS, double Ef, double int_Emin)
            : _pT(pT), _pS(pS), _Ef(Ef), _int_Emin(int_Emin) {}
        double operator() () const;
    private:
        const Sptr_Transmission_probability _pT;
        const Sptr_Electron_supply _pS;
        double _Ef;
        double _int_Emin;
};

#endif // _CURRENT_H_
