#include "../../include/distributors/Plummer.h"

// http://articles.adsabs.harvard.edu/pdf/1974A%26A....37..183A
Plummer::Plummer(unsigned long seed, int numParticles) : Distributor(seed), numParticles(numParticles) {

    confP = ConfigParser("config/Plummer.info");
    std::string description = confP.getVal<std::string>("description");

    std::cout << "Description: " << description << std::endl;

    M = confP.getVal<double>("M");
    R = confP.getVal<double>("R");
    R_max = confP.getVal<double>("R_max");

    rnd1 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd2 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd3 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd4 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd5 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd6 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd7 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));

}

Particle Plummer::next(int i) {

    Particle particle;
    double r;
    do {
         r = pow(pow(rnd1(gen), -2./3.) - 1., -1./2.);
    } while (r*R > R_max);

    double z = (1. - 2.*rnd2(gen)) * r;
    double rnd3_temp = rnd3(gen);
    double x = sqrt(r*r - z*z) * cos(2. * M_PI * rnd3_temp);
    double y = sqrt(r*r - z*z) * sin(2. * M_PI * rnd3_temp);

    particle.pos = vec3(x ,y, z) * R;

    double v_e = sqrt(2.) * pow(1. + r*r, -1./4.); // escape velocity

    double rnd4_temp;
    double rnd5_temp;

    do {
        rnd4_temp = rnd4(gen);
        rnd5_temp = rnd5(gen);
    } while (0.1 * rnd5_temp >= rnd4_temp*rnd4_temp * pow(1.-rnd4_temp*rnd4_temp, 7./2.));

    double v_mag = rnd4_temp * v_e;

    double vz = (1. - 2.*rnd6(gen)) * v_mag;
    double rnd7_temp = rnd7(gen);
    double vx = sqrt(v_mag * v_mag - vz * vz) * cos(2. * M_PI * rnd7_temp);
    double vy = sqrt(v_mag * v_mag - vz * vz) * sin(2. * M_PI * rnd7_temp);

    particle.vel = vec3(vx, vy, vz)  * sqrt(M/R);

    particle.mass = M/numParticles;

    return particle;
}
