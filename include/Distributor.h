#ifndef PARTICLEDISTRIBUTOR_DISTRIBUTOR_H
#define PARTICLEDISTRIBUTOR_DISTRIBUTOR_H

//#include "distributors/Plummer.h"
#include "Particle.h"

#include <random>
#include <iostream>

//#define _USE_MATH_DEFINES
#include <cmath>

class Distributor {

private:

    unsigned long seed;

protected:

    std::mt19937 gen;
    //const double G = 6.67408e-11;
    const double G = 1.;

public:

    Distributor(unsigned long seed = 0UL);
    virtual ~Distributor() {};

    virtual Particle next() = 0;
    //static std::string getName() { return std::string{ "distributor" }; };
};


#endif //PARTICLEDISTRIBUTOR_DISTRIBUTOR_H
