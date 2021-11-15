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

    plummerCenterVec = new vec3[numPlummerSpheres];
    for (int n=0; n<numPlummerSpheres; ++n){
        // Plummer spheres centers shall have the same locations for the same random seed
        // and therefore centers have to be created on initialization to ensure this
        // behaviour for different numbers of particles as gen is used for all random distributions
        plummerCenterVec[n] = vec3(rndX(gen), rndY(gen), rndZ(gen));
    }

    currentPlummerCenter = plummerCenterVec;
    std::cout << "Creating Plummer sphere @" << *currentPlummerCenter << std::endl;
    particleCounter = 0;
}

MultiplePlummer::~MultiplePlummer(){
    delete[] plummerCenterVec;
}

Particle MultiplePlummer::next(int i){

    if (particleCounter > 0 && particleCounter % (numParticles/numPlummerSpheres) == 0){
        ++currentPlummerCenter;
        std::cout << "Creating Plummer sphere @" << *currentPlummerCenter << std::endl;
    }
    Particle p_ = Plummer::next(i);
    p_.pos += *currentPlummerCenter;
    p_.mass = p_.mass * numPlummerSpheres; // m=M/numParticles in base class
    ++particleCounter;
    return p_;
}