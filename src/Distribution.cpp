#include "../include/Distribution.h"

Distribution::Distribution(int numParticles, DistributionType distributionType, unsigned long seed) :
                numParticles(numParticles), seed(seed) {

    if (this->seed == 0UL) {
        std::random_device rd; // obtain a random number from hardware
        this->seed = rd();
    }

    switch (distributionType) {
        case DistributionType::plummer: {
            distributor = new Plummer(seed, numParticles);
        } break;
        case DistributionType::singleSpiral: {
            distributor = new SingleSpiral(seed, numParticles);
        } break;
        case DistributionType::kepler: {
            distributor = new Kepler(seed, numParticles);
        } break;
        case DistributionType::multiplePlummer: {
            distributor = new MultiplePlummer(seed, numParticles);
        } break;
        default:
            printf("not implemented!\n");
            exit(0);
    }

    std::cout << "Selected distribution type: " << distributor->getName() << std::endl;

}

Distribution::~Distribution() {

    delete distributor;

}

void Distribution::generate() {
    for (int i=0; i<numParticles; i++) {
        particles.push_back(distributor->next(i));
    }
}

void Distribution::write2file(const std::string& filename) {

    HighFive::File file("output/" + filename + "N" + std::to_string(numParticles) + "seed" + std::to_string(seed) + ".h5",
              HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    std::vector<double> _mass;
    std::vector<std::vector<double>> _pos;
    std::vector<std::vector<double>> _vel;

    for (int i=0; i<numParticles; i++) {
        _mass.push_back(particles[i].mass);
        _pos.push_back({ particles[i].pos.x, particles[i].pos.y, particles[i].pos.z });
        _vel.push_back({ particles[i].vel.x, particles[i].vel.y, particles[i].vel.z });
    }

    // create data sets
    HighFive::DataSet mass = file.createDataSet<double>("/m", HighFive::DataSpace::From(_mass));
    HighFive::DataSet pos = file.createDataSet<double>("/x",  HighFive::DataSpace::From(_pos));
    HighFive::DataSet vel = file.createDataSet<double>("/v",  HighFive::DataSpace::From(_vel));


    // write data
    mass.write(_mass);
    pos.write(_pos);
    vel.write(_vel);

}