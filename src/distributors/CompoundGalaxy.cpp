#include "../../include/distributors/CompoundGalaxy.h"

CompoundGalaxy::CompoundGalaxy(unsigned long seed, int numParticles) : Distributor(seed), numParticles(numParticles) {
    confP = ConfigParser("config/CompoundGalaxy.info");
    std::string description = confP.getVal<std::string>("description");

    std::cout << "Description: " << description << std::endl;

    // disk parameters
    M_d = confP.getVal<double>("M_d");
    h = confP.getVal<double>("h");
    z_0 = confP.getVal<double>("z_0");

    // bulge parameters
    M_b = confP.getVal<double>("M_b");
    a = confP.getVal<double>("a");
    c = confP.getVal<double>("c");

    // halo parameters
    M_h = confP.getVal<double>("M_h");
    r_c = confP.getVal<double>("r_c");
    gamma = confP.getVal<double>("gamma");
    q = gamma/r_c;
    alpha = 1./(1- sqrt(M_PI)*q*exp(q*q)*(1-erf(q)));
    std::cout << "    Normalization constant for halo: alpha = " << alpha << std::endl;

    rnd1 = std::uniform_real_distribution<double>(0., std::nextafter(M_d, std::numeric_limits<double>::max()));
    rnd2 = std::uniform_real_distribution<double>(0., std::nextafter(M_d, std::numeric_limits<double>::max()));
    rnd3 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd4 = std::uniform_real_distribution<double>(0., std::nextafter(M_b, std::numeric_limits<double>::max()));
    rnd5 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd6 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd7 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd8 = std::uniform_real_distribution<double>(0., std::nextafter(M_h, std::numeric_limits<double>::max()));
    rnd9 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd10 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd11 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd12 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd13 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd14 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd15 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));

    // equal number of particles for each component
    // TODO: fix to work for arbitrary particle counts maybe by additional input parameters
    if (numParticles % 3 != 0){
        std::cerr << "Please select N as a multiple of three. - Aborting." << std::endl;
        exit(0);
    }
    // create positions of all components first for mass integration i velocity computation
    diskParticles = std::vector<Particle>(numParticles/3);
    bulgeParticles = std::vector<Particle>(numParticles/3);
    haloParticles = std::vector<Particle> (numParticles/3);
    // store radius in vector to achieve mass integration by index finding
    // equal mass particles per component is necessary for this method
    diskRadii = std::vector<double>(numParticles/3);
    bulgeRadii = std::vector<double>(numParticles/3);
    haloRadii = std::vector<double>(numParticles/3);

    for(int i=0; i<numParticles/3;++i){
        posDisk(diskParticles[i]);
        posBulge(bulgeParticles[i]);
        std::cout << "  Halo particle " << i << std::endl;
        posHalo(haloParticles[i]);
        diskRadii[i] = diskParticles[i].pos.getMagnitude();
        bulgeRadii[i] = bulgeParticles[i].pos.getMagnitude();
        haloRadii[i] = haloParticles[i].pos.getMagnitude();
    }

    // sort radii vectors
    std::sort(diskRadii.begin(), diskRadii.end());
    std::sort(bulgeRadii.begin(), bulgeRadii.end());
    std::sort(haloRadii.begin(), haloRadii.end());

}

Particle CompoundGalaxy::next(int i) {
    Particle particle;
    switch (component) {
        case disk: {
            particle = diskParticles[i];
            velDisk(particle);
        } break;
        case bulge: {
            particle = bulgeParticles[i-numParticles/3];
            velBulge(particle);
        } break;
        case halo: {
            particle = haloParticles[i-2*numParticles/3];
            std::cout << "Sampling velocity of particle " << i << std::endl;
            velHalo(particle);
        } break;
        default:
            std::cerr << "Unknown component type encountered. - Aborting." << std::endl;
            exit(0);
    }
    ++particleCounter;
    if (particleCounter % (numParticles/3) == 0){
        ++component;
    }
    return particle;
}

void CompoundGalaxy::posDisk(Particle &p) {
    //double M = rnd1(gen);
    double R = - h * log(1 - h/M_d*rnd1(gen));
    double z = z_0 * atanh(2*h/M_d*rnd2(gen)-1.);

    double phi = 2. * M_PI * rnd3(gen);
    double x = R * cos(phi);
    double y = R * sin(phi);

    p.pos = vec3(x, y, z);
    p.mass = M_d/(numParticles/3);
    p.matId = disk;
}

void CompoundGalaxy::posBulge(Particle &p) {

    double x_, y_, z_;

    do {
        // draw random point on sphere with radius 1
        z_ = (1. - 2.*rnd5(gen));
        double phi = 2. * M_PI * rnd6(gen);
        x_ = sqrt(1. - z_*z_) * cos(phi);
        y_ = sqrt(1. - z_*z_) * sin(phi);

        // coordinate transformation and rejection of points
        // create a uniform distribution on the selected ellipsoidal surface
        // https://math.stackexchange.com/a/982833

    } while(rnd7(gen) >= sqrt(pow(a*c*y_, 2.) + pow(a*a*z_, 2.) + pow(a*c*x_, 2.))/(a*a));

    // elliptical radius m
    double m = 1./(sqrt(M_b/rnd4(gen)) - 1.);

    p.pos = vec3(a*x_, a*y_, c*z_)*m;
    p.mass = M_b/(numParticles/3);
    p.matId = bulge;
}

void CompoundGalaxy::posHalo(Particle &p) {

    double M = rnd8(gen);

    auto f4root = [&](double _r){
        double rho = M_h*alpha*exp(-_r*_r/(r_c*r_c))/(2.*pow(M_PI, 3./2.)*r_c*(_r*_r+gamma*gamma));

        auto f2iMass = [&](double x){ return x*x*exp(-x*x)/(x*x+q*q); };
        double mass = 2*M_h*alpha/sqrt(M_PI)*boost::math::quadrature::trapezoidal(f2iMass, 0., _r/r_c);

        return std::pair<double,double>(mass - M, rho);
    };

    // arbitrary factor of 15*r_c as maximum value for r
    double r = boost::math::tools::newton_raphson_iterate(f4root, M/M_h*r_c, 0., 15.*r_c, 15);

    double z = (1. - 2.*rnd10(gen)) * r;
    double phi = 2. * M_PI * rnd11(gen);
    double x = sqrt(r*r - z*z) * cos(phi);
    double y = sqrt(r*r - z*z) * sin(phi);

    p.pos = vec3(x, y, z);
    p.mass = M_h/(numParticles/3);
    p.matId = halo;
}

void CompoundGalaxy::velDisk(Particle &p){
    //TODO: Implement
}

void CompoundGalaxy::velBulge(Particle &p){
    //TODO: Implement
}

void CompoundGalaxy::velHalo(Particle &p){
    double r = p.pos.getMagnitude();

    auto f2iHalo = [&](double _r)
            { return M_h*alpha*exp(-_r*_r/(r_c*r_c))/(2.*pow(M_PI, 3./2.)*r_c*(_r*_r+gamma*gamma))*G*massInSphere(_r)/(_r*_r); };

    // mean velocity @r squared
    double v2_r = 2.*pow(M_PI, 3./2.)*r_c*(r*r + gamma*gamma)/(M_h*alpha*exp(-r*r/(r_c*r_c)))
                    * boost::math::quadrature::trapezoidal(f2iHalo, r, 10.*r_c);

    //double sigma2 = pow(v2_r, 9./4.)/(8.*pow(2., 1./4.)*pow(M_PI, 13./4.)); // normalization factor squared
    double sigma2 = v2_r;

    // choose F between 0 and 1
    double F = rnd12(gen);// * 4.*M_PI*pow(2.*M_PI*sigma2, -3./2.)*2.*v2_r/M_E;

    std::cout << "  factor=" << 4.*M_PI*pow(2.*M_PI*sigma2, -3./2.) << std::endl;
    std::cout << "  r=" << r << ", M(r)=" << massInSphere(r) << ", F=" << F << ", sigma2=" << sigma2 << ", v2_r=" << v2_r << std::endl;

    auto f4root = [&](double _v){

        auto f2iF = [&](double x){ return x*x*exp(-x*x/(2.*v2_r)); };

        double F_ = 4.*M_PI*pow(2.*M_PI*sigma2, -3./2.)*boost::math::quadrature::trapezoidal(f2iF, 0., _v);
        double dFdv = 4.*M_PI*pow(2.*M_PI*sigma2, -3./2.)*_v*_v*exp(-_v*_v/(2.*v2_r));

        return std::pair<double, double>(F_ - F, dFdv);
    };

    // arbitrarily limit speed to 15*v_esc
    double v = boost::math::tools::newton_raphson_iterate(f4root, sqrt(v2_r), 0., 5.*sqrt(2.*G*massInSphere(r)/r), 15);

    std::cout << "Chosen v = " << v << std::endl;
    //v = 1./(4.*M_PI)*pow(2.*M_PI*sigma2, 3./2.)*v; // scale v

    // isotropy
    double v_z = (1. - 2.*rnd13(gen)) * v;
    double phi = 2. * M_PI * rnd14(gen);
    double v_x = sqrt(v*v - v_z*v_z) * cos(phi);
    double v_y = sqrt(v*v - v_z*v_z) * sin(phi);

    p.vel = vec3(v_x, v_y, v_z);
}

double CompoundGalaxy::massInSphere(double r){
    // greater than r predicate
    auto gtr = [&r](double val){ return val > r; };

    // approximating mass contained by finding index of r in sorted radii vectors
    // disk mass
    auto ir_d = std::find_if(begin(diskRadii), end(diskRadii), gtr);
    // bulge mass
    auto ir_b = std::find_if(begin(bulgeRadii), end(bulgeRadii), gtr);
    // halo mass
    auto ir_h = std::find_if(begin(haloRadii), end(haloRadii), gtr);

    return (double)(ir_d - diskRadii.begin())/(double)(numParticles/3)*M_d
            + (double)(ir_b - bulgeRadii.begin())/(double)(numParticles/3)*M_b
            + (double)(ir_h - haloRadii.begin())/(double)(numParticles/3)*M_h;
}