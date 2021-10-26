#ifndef PARTICLEDISTRIBUTOR_PLUMMER_H
#define PARTICLEDISTRIBUTOR_PLUMMER_H

#include "../Distributor.h"

class Plummer : public Distributor {

public:
    Plummer(unsigned long seed, int numParticles);

    Particle next(int i=0) override;

    static constexpr const char* name = "plummer";
    const std::string& getName() const override {
        static std::string name_ { name };
        return name_;
    }

protected:
    int numParticles;

    // values read from config file
    double M;

private:
    double R;
    double R_max;

    std::uniform_real_distribution<double> rnd1;
    std::uniform_real_distribution<double> rnd2;
    std::uniform_real_distribution<double> rnd3;
    std::uniform_real_distribution<double> rnd4;
    std::uniform_real_distribution<double> rnd5;
    std::uniform_real_distribution<double> rnd6;
    std::uniform_real_distribution<double> rnd7;
};


#endif //PARTICLEDISTRIBUTOR_PLUMMER_H
