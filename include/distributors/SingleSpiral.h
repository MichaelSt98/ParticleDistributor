#ifndef PARTICLEDISTRIBUTOR_SINGLESPIRAL_H
#define PARTICLEDISTRIBUTOR_SINGLESPIRAL_H

#include "../Distributor.h"

class SingleSpiral : public Distributor {

private:

    std::uniform_real_distribution<double> rndRCube; //(0., std::nextafter(1., std::numeric_limits<double>::max()));
    std::uniform_real_distribution<double> rndPhi; //(0., 2.*M_PI);
    std::uniform_real_distribution<double> rndCosTheta; //(-1., std::nextafter(1., std::numeric_limits<double>::max()));

public:

    SingleSpiral(unsigned long seed = 0UL);

    Particle next();
    static std::string getName() { return std::string{ "single_spiral" }; };
};


#endif //PARTICLEDISTRIBUTOR_SINGLESPIRAL_H
