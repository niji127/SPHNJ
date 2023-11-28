#ifndef _PARAMETER_H
#define _PARAMETER_H

#include "config.h"
#include "global.h"
#include "properties.h"
#include "physics.h"
#include "domain.h"
#include "sim_time.h"
#include "kernel.h"

#include "Eigen/Dense"

#include <string>
#include <fstream>

class Parameter
{
public:
    Physics physics;
    FluidProperties fluid;
    SolidProperties solid;
    VirtualProperties virt;
    Domain domain;
    SimTime time;
    Kernel kernel;

    void initiate()
    {
        readParticleNumber();
        readFieldParameter();
        initiateParameter();
    }

    void showInfo();
    void initiateParameter();

    // read parameters in input.ini
    void readFieldParameter();

    // read particle number (different files)
    void readParticleNumber();

    bool hasFluid();
    bool hasSolid();
    bool hasVirtual();

    SimTime *getTime();
};

#endif