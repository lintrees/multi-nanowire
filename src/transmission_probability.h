#ifndef _TRANSMISSION_PROBABILITY_H_
#define _TRANSMISSION_PROBABILITY_H_

#include "potential.h"

class Transmission_probability
{
    public:
        Transmission_probability(const Potential_energy U): _pU(&U) {};
        double operator() (double E) const;
    private:
        const Potential_energy* _pU;
};


#endif // _TRANSMISSION_PROBABILITY_H_
