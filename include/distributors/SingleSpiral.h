#ifndef PARTICLEDISTRIBUTOR_SINGLESPIRAL_H
#define PARTICLEDISTRIBUTOR_SINGLESPIRAL_H

#include "../Distributor.h"

class SingleSpiral : public Distributor {

private:

    int numParticles;
    std::uniform_real_distribution<double> rndRCube; //(0., std::nextafter(1., std::numeric_limits<double>::max()));
    std::uniform_real_distribution<double> rndPhi; //(0., 2.*M_PI);
    std::uniform_real_distribution<double> rndCosTheta; //(-1., std::nextafter(1., std::numeric_limits<double>::max()));

public:

    SingleSpiral(unsigned long seed, int numParticles);

    Particle next(int i=0) override;

    static constexpr const char* name = "single_spiral";

    const std::string& getName() const override {
        static std::string name_ { name };
        return name_;
    }

private:
    double M; // total mass
    double R; // radius of sphere
};


#endif //PARTICLEDISTRIBUTOR_SINGLESPIRAL_H
