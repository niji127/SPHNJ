#include "calculation.cuh"
#include "particle.cuh"
#include "memory_copy.cuh"

#include "thrust/iterator/constant_iterator.h"
#include "thrust/host_vector.h"

#include <fstream>

void Fluid_List::sort()
{
    auto fluid_prop_zip = getZipProperties();
    thrust::sort_by_key(thrust::device, cell_id.begin(), cell_id.end(), fluid_prop_zip);

    cell_id_list = cell_id;
    auto new_end = thrust::reduce_by_key(thrust::device, cell_id_list.begin(), cell_id_list.end(), thrust::constant_iterator<int>(1), cell_id_list.begin(), cell_counts.begin());

    cell_num = new_end.first - cell_id_list.begin();
    thrust::exclusive_scan(cell_counts.begin(), cell_counts.begin() + cell_num, cell_indices.begin());
}

void Solid_List::sort()
{
    auto solid_prop_zip = getZipProperties();
    cell_id_list = cell_id;
    thrust::sort_by_key(thrust::device, cell_id.begin(), cell_id.end(), solid_prop_zip);
    thrust::sort_by_key(thrust::device, cell_id_list.begin(), cell_id_list.end(), max_strain.begin());

    auto new_end = thrust::reduce_by_key(thrust::device, cell_id_list.begin(), cell_id_list.end(), thrust::constant_iterator<int>(1), cell_id_list.begin(), cell_counts.begin());

    cell_num = new_end.first - cell_id_list.begin();
    thrust::exclusive_scan(cell_counts.begin(), cell_counts.begin() + cell_num, cell_indices.begin());
}

void Virtual_List::sort()
{
    auto virt_prop_zip = getZipProperties();
    thrust::sort_by_key(thrust::device, cell_id.begin(), cell_id.end(), virt_prop_zip);

    cell_id_list = cell_id;
    auto new_end = thrust::reduce_by_key(thrust::device, cell_id_list.begin(), cell_id_list.end(), thrust::constant_iterator<int>(1), cell_id_list.begin(), cell_counts.begin());

    cell_num = new_end.first - cell_id_list.begin();
    thrust::exclusive_scan(cell_counts.begin(), cell_counts.begin() + cell_num, cell_indices.begin());
}

__global__ void clearCellNumber()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.domain.cell_number_total)
        return;

    Cell *cell = dev_data.cell + id;
    cell->fluid_num = 0;
    cell->solid_num = 0;
    cell->virt_num = 0;
}

__global__ void updateCellFluid(int *cell_id, int *cell_counts, int *cell_indices, const int fluid_cell_num)
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= fluid_cell_num)
        return;
    Cell *cell = dev_data.cell + cell_id[id];
    cell->fluid_num = cell_counts[id];
    for (int p_id = 0; p_id < cell_counts[id]; p_id++)
        cell->particle_list[p_id] = cell_indices[id] + p_id;
}

__global__ void updateCellSolid(int *cell_id, int *cell_counts, int *cell_indices, const int solid_cell_num)
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= solid_cell_num)
        return;

    Cell *cell = dev_data.cell + cell_id[id];
    cell->solid_num = cell_counts[id];
    for (int p_id = 0; p_id < cell_counts[id]; p_id++)
        cell->particle_list[p_id + cell->fluid_num] = cell_indices[id] + p_id;
}

__global__ void updateCellVirt(int *cell_id, int *cell_counts, int *cell_indices, const int virt_cell_num)
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= virt_cell_num)
        return;

    Cell *cell = dev_data.cell + cell_id[id];
    cell->virt_num = cell_counts[id];
    for (int p_id = 0; p_id < cell_counts[id]; p_id++)
        cell->particle_list[p_id + cell->fluid_num + cell->solid_num] = cell_indices[id] + p_id;
}

void updateCell(Sim *sim)
{
    int thread_num = Sim::parameter.kernel.thread_num;
    int cell_num = Sim::parameter.domain.cell_number_total;
    int cell_block = cell_num / thread_num + 1;
    clearCellNumber<<<cell_block, thread_num>>>();
    cudaDeviceSynchronize();
    if (Sim::parameter.hasFluid())
    {
        int *cell_id = thrust::raw_pointer_cast(sim->fluid_list.cell_id_list.data());
        int *cell_counts = thrust::raw_pointer_cast(sim->fluid_list.cell_counts.data());
        int *cell_indices = thrust::raw_pointer_cast(sim->fluid_list.cell_indices.data());

        int fluid_cell_num = sim->fluid_list.cell_num;
        int fluid_cell_block = fluid_cell_num / thread_num + 1;

        updateCellFluid<<<fluid_cell_block, thread_num>>>(cell_id, cell_counts, cell_indices, fluid_cell_num);
        cudaDeviceSynchronize();
    }

    if (Sim::parameter.hasSolid())
    {
        int *cell_id = thrust::raw_pointer_cast(sim->solid_list.cell_id_list.data());
        int *cell_counts = thrust::raw_pointer_cast(sim->solid_list.cell_counts.data());
        int *cell_indices = thrust::raw_pointer_cast(sim->solid_list.cell_indices.data());

        int solid_cell_num = sim->solid_list.cell_num;
        int solid_cell_block = solid_cell_num / thread_num + 1;
        updateCellSolid<<<solid_cell_block, thread_num>>>(cell_id, cell_counts, cell_indices, solid_cell_num);
        cudaDeviceSynchronize();
    }

    if (Sim::parameter.hasVirtual())
    {
        int *cell_id = thrust::raw_pointer_cast(sim->virt_list.cell_id_list.data());
        int *cell_counts = thrust::raw_pointer_cast(sim->virt_list.cell_counts.data());
        int *cell_indices = thrust::raw_pointer_cast(sim->virt_list.cell_indices.data());

        int virt_cell_num = sim->virt_list.cell_num;
        int virt_cell_block = virt_cell_num / thread_num + 1;
        updateCellVirt<<<virt_cell_block, thread_num>>>(cell_id, cell_counts, cell_indices, virt_cell_num);
        cudaDeviceSynchronize();
    }
}

void DeviceCalculation::sortParticle(Sim *sim)
{
    SimTime *time = Sim::parameter.getTime();
    if (!(time->isSort()))
        return;
    if (Sim::parameter.hasFluid())
        sim->fluid_list.sort();
    if (Sim::parameter.hasSolid())
        sim->solid_list.sort();
    if (Sim::parameter.hasVirtual())
        sim->virt_list.sort();
    updateCell(sim);
}