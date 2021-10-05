#ifndef PARTICLEDISTRIBUTOR_PLUMMER_H
#define PARTICLEDISTRIBUTOR_PLUMMER_H

#include "../Distributor.h"

class Plummer : public Distributor {

private:

    std::uniform_real_distribution<double> rnd1;
    std::uniform_real_distribution<double> rnd2;
    std::uniform_real_distribution<double> rnd3;
    std::uniform_real_distribution<double> rnd4;
    std::uniform_real_distribution<double> rnd5;
    std::uniform_real_distribution<double> rnd6;
    std::uniform_real_distribution<double> rnd7;

    int numParticles;

public:

    Plummer(unsigned long seed = 0UL, int numParticles=100000);

    Particle next();
    static std::string getName() { return std::string{ "plummer" }; };

};


#endif //PARTICLEDISTRIBUTOR_PLUMMER_H
