#include "../../include/distributors/Kepler.h"

Kepler::Kepler(unsigned long seed, int numParticles) : Distributor(seed){

    confP = ConfigParser("config/Kepler.info");
    std::string description = confP.getVal<std::string>("description");
    double diskWidth = confP.getVal<double>("diskWidth");

    std::cout << "Description: " << description << std::endl;

    R = confP.getVal<double>("R");
    M = confP.getVal<double>("M");
    m = confP.getVal<double>("m");
    applyPerturbations = confP.getVal<bool>("perturbations");
    perturbationDim = confP.getVal<double>("perturbationDim");

    perturbations = std::uniform_int_distribution<unsigned long>(0, 1000);
    distribution = std::uniform_real_distribution<double>(R - diskWidth, std::nextafter(R, std::numeric_limits<double>::max()));
    thetaDistribution = std::uniform_real_distribution<double>(0.0, 2 * M_PI);

}

Particle Kepler::next(int i) {

    Particle particle;

    double theta = thetaDistribution(gen);
    double r = distribution(gen);

    double G = 1; //6.67408e-11;

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
        if (applyPerturbations) {
            particle.pos = vec3{r * cos(theta), r * sin(theta),
                                (double) (perturbations(gen)) / 1000 * perturbationDim};
        }
        else {
            particle.pos = vec3{r * cos(theta), r * sin(theta), 0.0};
        }
        // velocity
        double rotation = 1;  // 1: clockwise   -1: counter-clockwise
        double v = sqrt(G * M / r);
        particle.vel = vec3{rotation * v * sin(theta), -rotation * v * cos(theta), 0.};
    }

    return particle;
}