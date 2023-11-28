#ifndef _OUTPUT_CUH
#define _OUTPUT_CUH

#include "function.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class Output : public DeviceFunction
{
public:
    static void writeFluidVelocity(Sim *sim, const std::string &path, SimTime *time);
    static void writeFluidDensity(Sim *sim, const std::string &path, SimTime *time);
    static void writeFluidPressure(Sim *sim, const std::string &path, SimTime *time);

    static void writeSolidVelocity(Sim *sim, const std::string &path, SimTime *time);
    static void writeSolidPressure(Sim *sim, const std::string &path, SimTime *time);
    static void writeSolidDensity(Sim *sim, const std::string &path, SimTime *time);
    static void writeSolidStress(Sim *sim, const std::string &path, SimTime *time);

    static void writeVirtNormal(Sim *sim, const std::string &path, SimTime *time);
    static void writeVirtPressure(Sim *sim, const std::string &path, SimTime *time);

    static void outputData(Sim *sim);

    static void outputEdge(Sim *sim, const std::string &path);
    static void outputTotalForce(Sim *sim, const std::string &path);
    static void outputWedgeDisp(Sim *sim, const std::string &path);
    static void outputMiddlePosition(Sim *sim, const std::string &path);
    static void outputEndPosition(Sim *sim, const std::string &path);
    static void outputTankBeamEnd(Sim *sim, const std::string &path);
    static void outputTotalEnergy(Sim *sim, const std::string &path);
    static void outputCrushForce(Sim *sim, const std::string &path);
    static void outputVonMisesStress(Sim *sim, const std::string &path);
    static void outputMeanPressure(Sim *sim, const std::string &path);
    static void outputGlobalData(Sim *sim);
};

#endif