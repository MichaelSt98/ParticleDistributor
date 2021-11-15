//
// Created by Johannes Martin on 21.10.21.
//

#ifndef PARTICLEDISTRIBUTOR_MULTIPLEPLUMMER_H
#define PARTICLEDISTRIBUTOR_MULTIPLEPLUMMER_H

#include "Plummer.h"

class MultiplePlummer : public Plummer {

public:

    MultiplePlummer(unsigned long seed, int numParticles);
    ~MultiplePlummer();

    Particle next(int i=0) override;

    static constexpr const char* name = "multiple_plummer";
    const std::string& getName() const override {
        static std::string name_ { name };
        return name_;
    }

private:

    std::uniform_real_distribution<double> rndX;
    std::uniform_real_distribution<double> rndY;
    std::uniform_real_distribution<double> rndZ;

    vec3 *plummerCenterVec;
    vec3 *currentPlummerCenter;

    // values read from config file
    int numPlummerSpheres;
    int particleCounter;
};


#endif //PARTICLEDISTRIBUTOR_MULTIPLEPLUMMER_H
