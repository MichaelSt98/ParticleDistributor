#include "../include/Distributor.h"

Distributor::Distributor(unsigned long seed){

    gen = std::mt19937(seed);
}