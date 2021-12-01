#ifndef PARTICLEDISTRIBUTOR_KEPLER_H
#define PARTICLEDISTRIBUTOR_KEPLER_H

#include "../Distributor.h"

using u32    = uint_least32_t;
using engine = std::mt19937;

class Kepler : public Distributor {

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
    bool applyPerturbations;
    double perturbationDim;

    std::uniform_int_distribution<unsigned long> perturbations;
    std::uniform_real_distribution<double> distribution;
    std::uniform_real_distribution<double> thetaDistribution;

};


#endif //PARTICLEDISTRIBUTOR_KEPLER_H
