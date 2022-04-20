#ifndef PARTICLEDISTRIBUTOR_COMPOUNDGALAXY_H
#define PARTICLEDISTRIBUTOR_COMPOUNDGALAXY_H

#include <algorithm>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/tools/roots.hpp>

#include "../Distributor.h"

class CompoundGalaxy : public Distributor {

public:
    CompoundGalaxy(unsigned long seed, int numParticles);

    Particle next(int i=0) override;

    static constexpr const char* name = "compound_galaxy";
    const std::string& getName() const override {
        static std::string name_ { name };
        return name_;
    }

protected:
    int numParticles;

    // values read from config file

    // disk related
    double M_d;
    double R_d;
    double h;
    double z_0;

    // bulge related
    double M_b;
    double a;
    double c;

    // halo related
    double M_h;
    double r_c;
    double gamma;
    // helper variables
    double q;
    double alpha;

private:

    enum Component { disk, bulge, halo };
    int component { disk };

    std::vector<Particle> diskParticles, bulgeParticles, haloParticles;
    std::vector<double> diskRadii, bulgeRadii, haloRadii;

    std::uniform_real_distribution<double> rnd1, rnd2, rnd3, rnd4, rnd5, rnd6, rnd7, rnd8, rnd9, rnd10, rnd11, rnd12,
                                           rnd13, rnd14, rnd15;

    int particleCounter { 0 };

    void posDisk(Particle &p);
    void posBulge(Particle &p);
    void posHalo(Particle &p);

    void velDisk(Particle &p);
    void velBulge(Particle &p);
    void velHalo(Particle &p);

    // functions to integrate numerically for velocity computation
    double massInSphere(double r);
};


#endif //PARTICLEDISTRIBUTOR_COMPOUNDGALAXY_H
