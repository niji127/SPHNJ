#ifndef _CUDA_MACRO_CUH
#define _CUDA_MACRO_CUH

#include <stdio.h>

// // disable stupid warnings
// #ifdef __INTELLISENSE__
// void __syncthreads();
// #define KERNEL_ARG2(grid, block)
// #define KERNEL_ARG3(grid, block, sh_mem)
// #define KERNEL_ARG4(grid, block, sh_mem, stream)
// #else
// #define KERNEL_ARG2(grid, block) <<< grid, block >>>
// #define KERNEL_ARG3(grid, block, sh_mem) <<< grid, block, sh_mem >>>
// #define KERNEL_ARG4(grid, block, sh_mem, stream) <<< grid, block, sh_mem, stream >>>
// #endif

// check cuda functions
#define CHECK(call)                                     \
    do                                                  \
    {                                                   \
        const cudaError_t error_code = call;            \
        if (error_code != cudaSuccess)                  \
        {                                               \
            printf("CUDA Error:\n");                    \
            printf("    File:       %s\n", __FILE__);   \
            printf("    Line:       %d\n", __LINE__);   \
            printf("    Error code: %d\n", error_code); \
            printf("    Error text: %s\n",              \
                   cudaGetErrorString(error_code));     \
            exit(1);                                    \
        }                                               \
    } while (0)

#define CUDA_ERROR_CHECK                             \
    if ((error = cudaGetLastError()) != cudaSuccess) \
        printf("CUDA error: %s\n", cudaGetErrorString(error));

// // overload for function atomicAdd (*double)
// #if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600

// #else
// __device__ double atomicAdd(double *address, double val)
// {
//     unsigned long long int *address_as_ull =
//         (unsigned long long int *)address;
//     unsigned long long int old = *address_as_ull, assumed;

//     do
//     {
//         assumed = old;
//         old = atomicCAS(address_as_ull, assumed,
//                         __double_as_longlong(val +
//                                              __longlong_as_double(assumed)));

//         // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
//     } while (assumed != old);

//     return __longlong_as_double(old);
// }
// #endif

#endif