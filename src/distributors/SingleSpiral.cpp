#include "../../include/distributors/SingleSpiral.h"

SingleSpiral::SingleSpiral(unsigned long seed, int numParticles) : Distributor(seed), numParticles(numParticles) {

    confP = ConfigParser("config/SingleSpiral.info");
    std::string description = confP.getVal<std::string>("description");

    std::cout << "Description: " << description << std::endl;

    rndRCube = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rndPhi = std::uniform_real_distribution<double>(0., 2.*M_PI);
    rndCosTheta =  std::uniform_real_distribution<double>(-1., std::nextafter(1., std::numeric_limits<double>::max()));

}

Particle SingleSpiral::next(int i) {

    Particle particle;
    double M = 1;
    double R = 1;

    double r, phi, theta, vmag;

    // generate uniformly distributed particles in a sphere
    r = R * pow(rndRCube(gen), 1./3.);
    phi = rndPhi(gen);
    theta = acos(rndCosTheta(gen));

    vmag = sqrt(G*M*r*r/(2.*R*R*R));

    particle.pos = vec3{ r * cos(phi) * sin(theta),
                         r * sin(phi) * sin(theta),
                         r * cos(theta) };

    particle.vel = vec3{ vmag * sin(phi),
                         -vmag * cos(phi),
                         0.};

    particle.mass = M/numParticles;

    //std::cout << "Particle.pos: " << particle.pos;

    return particle;

}