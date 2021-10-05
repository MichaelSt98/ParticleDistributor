#include "../include/Distribution.h"
#include "../include/Distributor.h"
#include "../include/Particle.h"

#include <cxxopts.hpp>
#include <sstream>

int main(int argc, char *argv[]) {

    cxxopts::Options options("ParticleDistributor",
                             "Generating HDF5 file with initial particle distribution for N-Body and SPH simulation.");

    std::vector<std::string> availableDistributionTypes = { Plummer::getName(), SingleSpiral::getName() };
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
            ("M,M-system", "Total mass distributed in the system", cxxopts::value<double>()->default_value("1."))
            ("R,R-sphere", "Radius of spheres serving as initial galaxies", cxxopts::value<double>()->default_value("1."))
            ("s,seed", "Use given random seed", cxxopts::value<unsigned long>()->default_value("0"))
            ("f,filename", "File name", cxxopts::value<std::string>()->default_value("test"))
            ("d,distributionType", availableDistributionTypesHelp.str(), cxxopts::value<int>()->default_value("0"))
            ("h,help", "Show this help");

    auto opts = options.parse(argc, argv);

    if (opts.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    std::string filename = opts["filename"].as<std::string>();

    const int numParticles { opts["N-particles"].as<int>() };
    const double M { opts["M-system"].as<double>() };
    const double R { opts["R-sphere"].as<double>() };

    DistributionType distributionType = DistributionType(opts["distributionType"].as<int>());
    unsigned long seed = opts["seed"].as<unsigned long>();

    Distribution distribution(numParticles, distributionType, seed);

    distribution.generate();
    distribution.write2file(filename);

    return 0;
}

