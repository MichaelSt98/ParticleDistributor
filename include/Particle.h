#ifndef PARTICLEDISTRIBUTOR_PARTICLE_H
#define PARTICLEDISTRIBUTOR_PARTICLE_H

#include "vector3.h"

typedef Vector3<double> vec3;

class Particle {

public:

    Particle();

    vec3 pos;
    vec3 vel;
    double mass;



};


#endif //PARTICLEDISTRIBUTOR_PARTICLE_H
