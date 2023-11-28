#include "sim.cuh"
#include "memory_copy.cuh"
#include "cuda_macro.cuh"

__constant__ Parameter dev_par;
__constant__ Data_Pointer dev_data;

void MemoryCopy::copyConstant()
{
    CHECK(cudaMemcpyToSymbol(dev_data, &(Sim::data_pointer), sizeof(Data_Pointer)));
    CHECK(cudaMemcpyToSymbol(dev_par, &(Sim::parameter), sizeof(Parameter)));
}