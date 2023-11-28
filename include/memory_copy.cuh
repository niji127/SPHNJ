#ifndef _MEMCOPY_CUH
#define _MEMCOPY_CUH

#include "cuda_macro.cuh"
#include "sim.cuh"
#include "function.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

extern __constant__ Parameter dev_par;
extern __constant__ Data_Pointer dev_data;

class MemoryCopy:public DeviceFunction
{
public:
    MemoryCopy() {}
    template <typename T>
    static T *hostToDevice(T *value, const int &num = 1, bool isDelete = false);

    template <typename T>
    static T **hostToDevice(T *hostData, T **devData, const int &rowNum, const int &colNum);

    static Sim *initiateGPU(Sim *host_sim);
    static void hostToDevice(Sim *host_sim, Sim *device_sim) {}
    static void deviceToHost(Sim *host_sim, Sim *device_sim) {}
    static void MemoryCopy::copyConstant();
};

template <typename T>
static T *MemoryCopy::hostToDevice(T *value, const int &num, bool isDelete)
{
    // A *m_a = new A[3];
    // A *a = valueHostToDevice(m_a, 3);
    T *devValue;
    CHECK(cudaMalloc((void **)&devValue, num * sizeof(T)));
    CHECK(cudaMemcpy(devValue, value, num * sizeof(T), cudaMemcpyHostToDevice));
    if (isDelete)
    {
        if (num == 1)
            delete value;
        else
            delete[] value;
    }
    return devValue;
}

template <typename T>
static T **MemoryCopy::hostToDevice(T *hostData, T **devData, const int &rowNum, const int &colNum)
{
    // vec **image2D = NULL, *gpu_image2DData = NULL;
    // vec *cpu_image2DData = NULL;
    // image2D = array2DHostToDevice(cpu_image2DData, &gpu_normalData, height, width);

    // allocate memory
    T **hostArray;
    cudaHostAlloc((void **)&hostArray, rowNum * sizeof(T *), cudaHostAllocDefault);
    T **devArray;
    cudaMalloc((void **)&devArray, rowNum * sizeof(T *));
    size_t pitch;
    cudaMallocPitch(devData, &pitch, colNum * sizeof(T), rowNum);
    cudaMemset2D(*devData, pitch, 0, colNum * sizeof(T), rowNum);
    // set pointer
    for (int i = 0; i < rowNum; i++)
        hostArray[i] = (T *)((float *)(*devData) + i * pitch / sizeof(float));
    // copy data
    if (hostData != NULL)
        cudaMemcpy2D(*devData, pitch, hostData, colNum * sizeof(T), colNum * sizeof(T), rowNum, cudaMemcpyHostToDevice);
    cudaMemcpy(devArray, hostArray, rowNum * sizeof(T *), cudaMemcpyHostToDevice);
    cudaFreeHost(hostArray);
    return devArray;
}

#endif