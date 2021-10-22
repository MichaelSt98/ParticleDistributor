#ifndef PARTICLEDISTRIBUTOR_KEPLER_H
#define PARTICLEDISTRIBUTOR_KEPLER_H

#include "../Distributor.h"

class Kepler : public Distributor {

    std::uniform_real_distribution<double> distribution;
    std::uniform_real_distribution<double> thetaDistribution;
    int numParticles;

public:

    Kepler(unsigned long seed, int numParticles);

    Particle next(int i=0) override;

    static constexpr const char* name = "kepler";
    const std::string& getName() const override {
        static std::string name_ { name };
        return name_;
    }

private:
    double R; // outer disk radius
    double M; // mass of central star
    double m; // mass of orbiting stars
};


#endif //PARTICLEDISTRIBUTOR_KEPLER_H
