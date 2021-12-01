#include "../../include/distributors/KeplerTorus.h"

KeplerTorus::KeplerTorus(unsigned long seed, int numParticles) : Distributor(seed){

    confP = ConfigParser("config/KeplerTorus.info");
    std::string description = confP.getVal<std::string>("description");

    std::cout << "Description: " << description << std::endl;

    R = confP.getVal<double>("R");
    r = confP.getVal<double>("r");
    M = confP.getVal<double>("M");
    m = confP.getVal<double>("m");

    phiDistribution = std::uniform_real_distribution<double>(0.0, 2 * M_PI);
    thetaDistribution = std::uniform_real_distribution<double>(0.0, 2 * M_PI);

}

Particle KeplerTorus::next(int i) {

    Particle particle;

    double theta = thetaDistribution(gen);
    double phi = phiDistribution(gen);

    double G = 1.; //6.67408e-11;

    if (i % 1000 == 0) { std::cout << "i = " << i << ": M = " << M << ", m = " << m << std::endl; }

    if (i == 0) {
        // mass
        particle.mass = M;
        // position
        particle.pos = vec3{ 0., 0., 0. };
        // velocity
        particle.vel = vec3{ 0., 0., 0. };
    }
    else {
        // mass
        particle.mass = m;
        // position
        particle.pos = vec3{(R + r * cos(theta)) * cos(phi),
                            (R + r * cos(theta)) * sin(phi),
                            r * sin(theta)};
        // velocity
        double rotation = 1;  // 1: clockwise   -1: counter-clockwise
        double _r = sqrt(pow(particle.pos.x, 2) +
                         pow(particle.pos.y, 2) +
                         pow(particle.pos.z, 2));

        double v = sqrt(G * M / _r);
        particle.vel = vec3{rotation * v * sin(phi), -rotation * v * cos(phi), 0.};
    }

    return particle;
}