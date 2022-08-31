#include "../include/Distribution.h"
#include "../include/Distributor.h"
#include "../include/Particle.h"

#include <cxxopts.hpp>
#include <sstream>

int main(int argc, char *argv[]) {

    cxxopts::Options options("bin/runner",
                             "Generating HDF5 file with initial particle distribution for N-Body and SPH simulation.");

    std::vector<std::string> availableDistributionTypes = { Plummer::name, SingleSpiral::name, Kepler::name,
	MultiplePlummer::name, KeplerTorus::name };
                                                            //,CompoundGalaxy::name };
							    // TODO: remove the above line when compound_galaxy distributor is working
    std::stringstream availableDistributionTypesHelp;

    availableDistributionTypesHelp << "(";
    for (int i=0; i<availableDistributionTypes.size(); i++) {
        if (i < availableDistributionTypes.size() - 1) {
            availableDistributionTypesHelp << "[" << i << "]: " << availableDistributionTypes[i] << ", ";
        }
        else {
            availableDistributionTypesHelp << "[" << i << "]: " << availableDistributionTypes[i] << ")";
        }
    }

    options.set_width(100);
    options.add_options()
            ("N,N-particles", "Number of particles", cxxopts::value<int>()->default_value("1000000"))
            ("s,seed", "Use given random seed", cxxopts::value<unsigned long>()->default_value("0"))
            ("f,filename", "File name", cxxopts::value<std::string>())
            ("d,distributionType", availableDistributionTypesHelp.str(), cxxopts::value<int>()->default_value("0"))
            ("h,help", "Show this help");

    auto opts = options.parse(argc, argv);

    if (opts.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    DistributionType distributionType = DistributionType(opts["distributionType"].as<int>());

    const int numParticles { opts["N-particles"].as<int>() };
    unsigned long seed = opts["seed"].as<unsigned long>();

    Distribution distribution(numParticles, distributionType, seed);

    distribution.generate();

    std::string filename;
    if (opts.count("filename")){
        filename = opts["filename"].as<std::string>();
    } else {
        filename = distribution.getDistributionType();
    }

    distribution.write2file(filename);

    return 0;
}

