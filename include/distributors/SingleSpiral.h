#ifndef PARTICLEDISTRIBUTOR_SINGLESPIRAL_H
#define PARTICLEDISTRIBUTOR_SINGLESPIRAL_H

#include "../Distributor.h"

class SingleSpiral : public Distributor {

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

    int numParticles;
    std::uniform_real_distribution<double> rndRCube;
    std::uniform_real_distribution<double> rndPhi;
    std::uniform_real_distribution<double> rndCosTheta;
};


#endif //PARTICLEDISTRIBUTOR_SINGLESPIRAL_H
