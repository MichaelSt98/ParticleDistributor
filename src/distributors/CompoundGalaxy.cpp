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

    rnd1 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd2 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd3 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd4 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd5 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd6 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd7 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
    rnd8 = std::uniform_real_distribution<double>(0., std::nextafter(1., std::numeric_limits<double>::max()));
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
    // helper variable m
    double m = 1./(sqrt(M_b/rnd4(gen)) - 1.);

    // random point on sphere with radius m
    double x_, y_, z_;

    do {
        // draw random point on sphere with radius m
        z_ = (1. - 2.*rnd5(gen)) * m;
        double phi = 2. * M_PI * rnd6(gen);
        x_ = sqrt(m*m - z_*z_) * cos(phi);
        y_ = sqrt(m*m - z_*z_) * sin(phi);

        // coordinate transformation and rejection of points
        // create a uniform distribution on the selected ellipsoidal surface
        // https://math.stackexchange.com/a/982833

    } while(rnd7(gen) >= sqrt(pow(a*c*y_, 2.) + pow(a*a*z_, 2.) + pow(a*c*x_, 2.))/(a*a));

    p.pos = vec3(a*x_, a*y_, c*z_);
    p.mass = M_b/(numParticles/3);
    p.matId = bulge;
}

void CompoundGalaxy::posHalo(Particle &p) {
    // rejection sampling of r
    double r;
    do {
        r = rnd8(gen);
    } while(rnd9(gen) >= gamma*gamma*exp(-r*r/(r_c*r_c))/(r*r+gamma*gamma));

    // scale r from [0, 1] to [0, r_c]
    r = r_c * r; // TODO: check if this is the right way to scale
    // r = 2.*pow(M_PI, 3./2.)*r_c/(M_h*alpha) * r;

    // draw random point on sphere with radius r
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

    auto f2iHalo = [&](double r_)
            { return M_h*alpha*exp(-r_*r_/(r_c*r_c))/(2.*pow(M_PI, 3./2.)*r_c*(r_*r_+gamma*gamma))*G*massInSphere(r_)/(r_*r_); };

    // mean velocity @r squared
    double v2_r = 2.*pow(M_PI, 3./2.)*r_c*(r*r + gamma*gamma)/(M_h*alpha*exp(-r*r/(r_c*r_c)))
                    * boost::math::quadrature::trapezoidal(f2iHalo, r, r_c);

    double v;
    do {
        // rejection sampling of v
        do {
            v = rnd12(gen);
        } while (rnd13(gen) >= M_E / (2. * v2_r) * v * v * exp(-v * v / (2. * v2_r * v2_r)));

        //double sigma2 = pow(v2_r, 9./4.)/(8.*pow(2., 1./4.)*pow(M_PI, 13./4.)); // normalization factor squared
        //v = 1./(4.*M_PI)*pow(2.*M_PI*sigma2, 3./2.)*v; // scale v
        v = sqrt(2.*G*massInSphere(r)/r) * v;
        std::cout << "    r = " << r << ", v = " << v << ", v_esc = " << sqrt(2.*G*massInSphere(r)/r) << ", v_r = " << sqrt(v2_r) << std::endl;
    } while(false); //v > .95 * sqrt(2.*G*massInSphere(r)/r)); // limit v to .95 of escape velocity

    // isotropy
    double v_z = (1. - 2.*rnd14(gen)) * v;
    double phi = 2. * M_PI * rnd15(gen);
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