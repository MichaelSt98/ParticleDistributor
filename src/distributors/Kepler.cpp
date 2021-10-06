#include "../../include/distributors/Kepler.h"

Kepler::Kepler(unsigned long seed, int numParticles, double width) : Distributor(seed), numParticles(numParticles) {

    distribution = std::uniform_real_distribution<double>(1. - width, std::nextafter(1., std::numeric_limits<double>::max()));
    thetaDistribution = std::uniform_real_distribution<double>(0.0, 2 * M_PI);

}

Particle Kepler::next(int i) {

    Particle particle;

    double M = 1.;
    double G = 1.;

    double mass = M/numParticles;
    double theta = thetaDistribution(gen);
    double r = distribution(gen);

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
        particle.mass = 0.01 * mass;
        // position
        particle.pos = vec3{r * cos(theta), r * sin(theta), 0.};
        // velocity
        double rotation = 1;  // 1: clockwise   -1: counter-clockwise
        double v = sqrt(M / r);
        particle.vel = vec3{rotation * v * sin(theta), -rotation * v * cos(theta), 0.};
    }

    return particle;
}