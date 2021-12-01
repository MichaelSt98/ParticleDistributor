#ifndef PARTICLEDISTRIBUTOR_DISTRIBUTION_H
#define PARTICLEDISTRIBUTOR_DISTRIBUTION_H

#include "Particle.h"

#include "distributors/Plummer.h"
#include "distributors/SingleSpiral.h"
#include "distributors/Kepler.h"
#include "distributors/MultiplePlummer.h"
#include "distributors/KeplerTorus.h"

#include <vector>
#include <string>
#include <highfive/H5File.hpp>
#include <iostream>

enum class DistributionType {
    plummer, singleSpiral, kepler, multiplePlummer, keplerTorus
};

class Distribution {

private:

    int numParticles;
    std::vector<Particle> particles;
    Distributor *distributor;
    unsigned long seed;

public:

    Distribution(int numParticles, DistributionType distributionType, unsigned long _seed=0UL);
    ~Distribution();

    void generate();
    void write2file(const std::string& filename);

    std::string getDistributionType() const {
        return distributor->getName();
    }
};


#endif //PARTICLEDISTRIBUTOR_DISTRIBUTION_H
