#include "calculation.cuh"
#include "memory_copy.cuh"
#include "global.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__device__ int findPositionCell(const Eigen::Vector2f &pos)
{
    float xmin = dev_par.domain.cell_xmin[0];
    float ymin = dev_par.domain.cell_xmin[1];
    float inter_inv_x = dev_par.domain.interval_inv[0];
    float inter_inv_y = dev_par.domain.interval_inv[1];
    int block_nx = dev_par.domain.cell_number[0];
    int block_ny = dev_par.domain.cell_number[1];

    return (int)((pos[1] - ymin) * inter_inv_y) + (int)((pos[0] - xmin) * inter_inv_x) * block_ny;
}

__global__ void fillFluidCellID()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int fluid_num = dev_par.fluid.number_total;
    if (id >= fluid_num)
        return;
    dev_data.fluid.cell_id[id] = findPositionCell(dev_data.fluid.position[id]);
}

__global__ void fillSolidCellID()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int solid_num = dev_par.solid.number;
    if (id >= solid_num)
        return;
    dev_data.solid.cell_id[id] = findPositionCell(dev_data.solid.position[id].cast<float>());
}

__global__ void fillVirtCellID()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int virt_num = dev_par.virt.number;
    if (id >= virt_num)
        return;
    dev_data.virt.cell_id[id] = findPositionCell(dev_data.virt.position[id]);
}

void DeviceCalculation::particleCellUpdate()
{
    int thread_num = Sim::parameter.kernel.thread_num;

    int fluid_num = Sim::parameter.fluid.number_total;
    int fluid_block = fluid_num / thread_num + 1;
    fillFluidCellID<<<fluid_block, thread_num>>>();

    int solid_num = Sim::parameter.solid.number;
    int solid_block = solid_num / thread_num + 1;
    fillSolidCellID<<<solid_block, thread_num>>>();

    int virt_num = Sim::parameter.virt.number;
    int virt_block = virt_num / thread_num + 1;
    fillVirtCellID<<<virt_block, thread_num>>>();

    cudaDeviceSynchronize();
}

__global__ void fillCellList()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.domain.cell_number_total)
        return;

    Cell *cell_i = dev_data.cell + id;
    int fluid_num(0), solid_num(0), virt_num(0);
    // int particle_list[CELL_SIZE];

    // fluid
    for (int n_id = 0; n_id < 9; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;
        int fluid_num_j = cell_j->fluid_num_previous;

        for (int p_id = 0; p_id < fluid_num_j; p_id++)
        {
            int fluid_id = cell_j->particle_list_previous[p_id];
            if (dev_data.fluid.cell_id[fluid_id] == id)
            {
                if (fluid_num >= CELL_SIZE)
                {
                    printf("too many particles in a cell fid=%d num=%d %d %d\n", id, fluid_num, solid_num, virt_num);
                    return;
                }
                cell_i->particle_list[fluid_num] = fluid_id;
                fluid_num++;
            }
        }
    }
    cell_i->fluid_num = fluid_num;

    // solid
    for (int n_id = 0; n_id < 9; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;
        int fluid_num_j = cell_j->fluid_num_previous;
        int solid_num_j = cell_j->solid_num_previous;
        for (int p_id = fluid_num_j; p_id < fluid_num_j + solid_num_j; p_id++)
        {
            int solid_id = cell_j->particle_list_previous[p_id];
            if (dev_data.solid.cell_id[solid_id] == id)
            {
                if (fluid_num + solid_num >= CELL_SIZE)
                {
                    printf("too many particles in a cell sid=%d num=%d %d %d\n", id, fluid_num, solid_num, virt_num);
                    return;
                }
                cell_i->particle_list[fluid_num + solid_num] = solid_id;
                solid_num++;
            }
        }
    }
    cell_i->solid_num = solid_num;

    // virt
    for (int n_id = 0; n_id < 9; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;
        int fluid_num_j = cell_j->fluid_num_previous;
        int solid_num_j = cell_j->solid_num_previous;
        int virt_num_j = cell_j->virt_num_previous;
        for (int p_id = fluid_num_j + solid_num_j; p_id < fluid_num_j + solid_num_j + virt_num_j; p_id++)
        {
            int virt_id = cell_j->particle_list_previous[p_id];
            if (dev_data.virt.cell_id[virt_id] == id)
            {
                if (fluid_num + solid_num + virt_num >= CELL_SIZE)
                {
                    printf("too many particles in a cell vid=%d num=%d %d %d\n", id, fluid_num, solid_num, virt_num);
                    return;
                }
                cell_i->particle_list[fluid_num + solid_num + virt_num] = virt_id;
                virt_num++;
            }
        }
    }
    cell_i->virt_num = virt_num;
}

__global__ void copyCellList()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.domain.cell_number_total)
        return;

    Cell *cell = dev_data.cell + id;
    cell->fluid_num_previous = cell->fluid_num;
    cell->solid_num_previous = cell->solid_num;
    cell->virt_num_previous = cell->virt_num;
    for (int i = 0; i < cell->fluid_num + cell->solid_num + cell->virt_num; i++)
        cell->particle_list_previous[i] = cell->particle_list[i];
}

void DeviceCalculation::fillCell()
{
    int thread_num = Sim::parameter.kernel.thread_num;
    int cell_num = Sim::parameter.domain.cell_number_total;
    int cell_block = cell_num / thread_num + 1;

    copyCellList<<<cell_block, thread_num>>>();
    cudaDeviceSynchronize();

    fillCellList<<<cell_block, thread_num>>>();
    cudaDeviceSynchronize();
}