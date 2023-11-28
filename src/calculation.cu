#include "function.cuh"
#include "calculation.cuh"
#include "memory_copy.cuh"
#include "sim.cuh"
#include "output.cuh"
#include "probe.h"

#include "thrust/host_vector.h"

#include <fstream>
#include <string>
#include <windows.h>

void DeviceFunction::GPUCalulation(Sim *sim)
{
    MemoryCopy::copyConstant();
    DeviceCalculation::calculation(sim);
}

void DeviceCalculation::calculation(Sim *sim)
{
    SimTime *time = Sim::parameter.getTime();

    clock_t start, end;

    // Probe *probe_0 = new Probe({5.366f, 0.19f}, 0);

    virtNormalCalculation();
    start = clock();
    while (time->addStep())
    {
        Output::outputData(sim);
        Output::outputGlobalData(sim);
        particleCellUpdate();
        sortParticle(sim);
        fillCell();
        coupleCalculation();
        virtUpdate();
        fluidCalculation();
        solidCalculation(sim);
        fluidCorrection(sim);
        // probe_0->outputData();
    }

    end = clock();
    printf("total time=%dms\n", end - start);
    Output::outputData(sim);
}