#ifndef _FUNCTION_CUH
#define _FUNCTION_CUH

#include "sim.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class DeviceFunction
{
public:
    static void GPUCalulation(Sim *sim);
    static void expandHeapSize();
    template <typename T>
    __device__ static void warpReduce(T *sdata, int tid);

    template <typename T>
    __device__ static void reduction(T *sdata, int tid, int block_size);
};

template <typename T>
__device__ void warpReduce(T *sdata, int tid)
{

    T v = sdata[tid];
    v += sdata[tid + 32];
    sdata[tid] = v;
    __syncwarp();
    v += sdata[tid + 16];
    __syncwarp();
    sdata[tid] = v;
    __syncwarp();
    v += sdata[tid + 8];
    __syncwarp();
    sdata[tid] = v;
    __syncwarp();
    v += sdata[tid + 4];
    __syncwarp();
    sdata[tid] = v;
    __syncwarp();
    v += sdata[tid + 2];
    __syncwarp();
    sdata[tid] = v;
    __syncwarp();
    v += sdata[tid + 1];
    __syncwarp();
    sdata[tid] = v;
}

template <typename T>
__device__ void reduction(T *sdata, int tid, int block_size)
{
    // uncomment to support blocksize of 1024
    if (block_size == 1024)
    {
        if (tid < 512)
            sdata[tid] += sdata[tid + 512];
        __syncthreads();
    }
    if (block_size >= 512)
    {
        if (tid < 256)
            sdata[tid] += sdata[tid + 256];
        __syncthreads();
    }
    if (block_size >= 256)
    {
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
    }
    if (block_size >= 128)
    {
        if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
    }
    if (tid < 32)
        warpReduce(sdata, tid);
}

#endif