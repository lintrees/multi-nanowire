#ifndef _CURRENT_H_
#define _CURRENT_H_

#include "transmission_probability.h"
#include "electron_supply.h"


class Current_density
{
    public:
        Current_density(const Transmission_probability& T,
            const Electron_supply& S, double Ef, double int_Emin)
            : _pT(&T), _pS(&S), _Ef(Ef), _int_Emin(int_Emin) {}
        double operator() () const;
    private:
        const Transmission_probability* _pT;
        const Electron_supply* _pS;
        double _Ef;
        double _int_Emin;
};

#endif // _CURRENT_H_
