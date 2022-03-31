#include "../include/Distribution.h"

Distribution::Distribution(int numParticles, DistributionType distributionType, unsigned long _seed) :
                numParticles(numParticles), distributionType(distributionType), seed(_seed) {

    if (seed == 0UL) {
        std::random_device rd; // obtain a random number from hardware
        seed = rd();
        std::cout << "Obtained random seed from hardware: " << seed << std::endl;
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
        case DistributionType::keplerTorus: {
            distributor = new KeplerTorus(seed, numParticles);
        } break;
        case DistributionType::compoundGalaxy: {
            distributor = new CompoundGalaxy(seed, numParticles);
        } break;
        default:
            std::cerr << "Not implemented. - Aborting." << std::endl;
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

    std::vector<double> _mass(numParticles);
    std::vector<std::vector<double>> _pos(numParticles);
    std::vector<std::vector<double>> _vel(numParticles);

    for (int i=0; i<numParticles; i++) {
        _mass[i] = particles[i].mass;
        _pos[i] = { particles[i].pos.x, particles[i].pos.y, particles[i].pos.z };
        _vel[i] = { particles[i].vel.x, particles[i].vel.y, particles[i].vel.z };
    }

    // create data sets
    HighFive::DataSet mass = file.createDataSet<double>("/m", HighFive::DataSpace::From(_mass));
    HighFive::DataSet pos = file.createDataSet<double>("/x",  HighFive::DataSpace::From(_pos));
    HighFive::DataSet vel = file.createDataSet<double>("/v",  HighFive::DataSpace::From(_vel));

    // adding material type if distribution is a compound galaxy
    if (distributionType == DistributionType::compoundGalaxy){
        std::vector<int> _matId(numParticles);
        for (int i=0; i<numParticles; i++) {
            _matId[i] = particles[i].matId;
        }
        HighFive::DataSet matId = file.createDataSet<int>("/materialId", HighFive::DataSpace::From(_matId));
        matId.write(_matId);
    }

    // write data
    mass.write(_mass);
    pos.write(_pos);
    vel.write(_vel);

}