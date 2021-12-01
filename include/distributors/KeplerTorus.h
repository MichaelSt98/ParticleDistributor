#ifndef PARTICLEDISTRIBUTOR_KEPLERTORUS_H
#define PARTICLEDISTRIBUTOR_KEPLERTORUS_H

#include "../Distributor.h"

class KeplerTorus : public Distributor {

public:
    KeplerTorus(unsigned long seed, int numParticles);

    Particle next(int i=0) override;

    static constexpr const char* name = "kepler_torus";
    const std::string& getName() const override {
        static std::string name_ { name };
        return name_;
    }

private:
    double r; // torus radius
    double R; // outer disk radius
    double M; // mass of central star
    double m; // mass of orbiting stars

    std::uniform_real_distribution<double> phiDistribution;
    std::uniform_real_distribution<double> thetaDistribution;

};


#endif //PARTICLEDISTRIBUTOR_KEPLER_H
