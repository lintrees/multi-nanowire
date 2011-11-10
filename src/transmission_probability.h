#ifndef _TRANSMISSION_PROBABILITY_H_
#define _TRANSMISSION_PROBABILITY_H_

#include <memory>

#include "potential.h"

namespace Cartesian_1d
{

class Transmission_probability
{
    public:
        Transmission_probability(const Sptr_Potential_energy pU): _pU(pU) {};
        double operator() (double E) const;
    private:
        const Sptr_Potential_energy _pU;
};
typedef std::shared_ptr<Transmission_probability> Sptr_Transmission_probability;

} // namespace Cartesian_1d

#endif // _TRANSMISSION_PROBABILITY_H_
