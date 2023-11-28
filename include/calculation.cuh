#ifndef _CALCULATION_CUH
#define _CALCULATION_CUH

#include "function.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class DeviceCalculation : public DeviceFunction
{
public:
    static void calculation(Sim *sim);

    static void particleCellUpdate();
    static void fillCell();
    static void sortParticle(Sim *sim);
    static void coupleCalculation();
    static void fluidCalculation();
    static void virtNormalCalculation();
    static void virtUpdate();
    static void solidCalculation(Sim *sim);
    static void fluidCorrection(Sim *sim);
};

#endif