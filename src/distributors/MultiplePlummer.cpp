//
// Created by Johannes Martin on 21.10.21.
//

#include "../../include/distributors/MultiplePlummer.h"

MultiplePlummer::MultiplePlummer(unsigned long seed, int numParticles) : Plummer(seed, numParticles) {

    confP = ConfigParser("config/MultiplePlummer.info");
    std::string description = confP.getVal<std::string>("description");

    std::cout << "Description: " << description << std::endl;

    numPlummerSpheres = confP.getVal<int>("numPlummerSpheres");
    // overwrite base class members
    M = confP.getVal<double>("M");
    R = confP.getVal<double>("R");
    R_max = confP.getVal<double>("R_max");

    if (numParticles % numPlummerSpheres != 0){
        std::cout << "WARNING: numPlummerSpheres =" << numPlummerSpheres
                  << " is not a divisor of numParticles N = " << numParticles;
        throw std::invalid_argument("numPlummerSpheres must be divisor of N!");
    }

    const double boxLength { confP.getVal<double>("boxLength") };
    rndX = std::uniform_real_distribution<double>(-.5*boxLength,
                                                 std::nextafter(.5*boxLength, std::numeric_limits<double>::max()));
    rndY = std::uniform_real_distribution<double>(-.5*boxLength,
                                                 std::nextafter(.5*boxLength, std::numeric_limits<double>::max()));
    rndZ = std::uniform_real_distribution<double>(-.5*boxLength,
                                                 std::nextafter(.5*boxLength, std::numeric_limits<double>::max()));

    particleCounter = 0;

}

Particle MultiplePlummer::next(int i){

    if (particleCounter % (numParticles/numPlummerSpheres) == 0){
        plummerCenterVec = vec3(rndX(gen), rndY(gen), rndZ(gen));
        std::cout << "Creating Plummer sphere @" << plummerCenterVec << std::endl;
    }
    Particle p_ = Plummer::next(i);
    p_.pos += plummerCenterVec;
    p_.mass = p_.mass * numPlummerSpheres; // m=M/numParticles in base class
    ++particleCounter;
    return p_;
}