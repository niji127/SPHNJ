#include "calculation.cuh"
#include "memory_copy.cuh"
#include "global.h"

extern __device__ float densityToPressure(float &density, const PhaseType &phase);
extern __device__ float pressureToDensity(float &pressure, const PhaseType &phase);

__global__ void getVirtNormal()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.virt.number)
        return;

    Virtual_Data virt = dev_data.virt;
    Solid_Data solid = dev_data.solid;
    Eigen::Vector2f x_i = virt.position[id];
    Eigen::Vector2f normal = Eigen::Vector2f::Zero();

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    for (int id_j = 0; id_j < dev_par.virt.number; id_j++)
    {
        Eigen::Vector2f dx = virt.position[id_j] - x_i;
        if (dx.squaredNorm() < impact_dis_2)
        {
            float q_2 = dx.squaredNorm() * h_inv_2;
            float kernel_coef = expf(-q_2) * coef * coef3to2;
            normal -= dx * kernel_coef;
        }
    }
    // for (int id_j = 0; id_j < dev_par.solid.number; id_j++)
    // {
    //     Eigen::Vector2f dx = solid.position[id_j].cast<float>() - x_i;
    //     if (dx.squaredNorm() < impact_dis_2)
    //     {
    //         float q_2 = dx.squaredNorm() * h_inv_2;
    //         float kernel_coef = expf(-q_2) * coef * coef3to2;
    //         normal -= dx * kernel_coef;
    //     }
    // }
    if (normal.norm() < EPS_FOR_SUM)
        virt.normal[id] = Eigen::Vector2f::Zero();
    else
        virt.normal[id] = normal.normalized();
}

void DeviceCalculation::virtNormalCalculation()
{
    if (!(Sim::parameter.hasVirtual()))
        return;
    int thread_num = Sim::parameter.kernel.thread_num;
    int virt_num = Sim::parameter.virt.number;
    int virt_block = virt_num / thread_num + 1;
    getVirtNormal<<<virt_block, thread_num>>>();
    cudaDeviceSynchronize();
}

__device__ void virtDataCalculate(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, Eigen::Vector2f &velocity, float &pressure, float &sum_correction, float &nearest_distance_square)
{
    Fluid_Data fluid = dev_data.fluid;
    Virtual_Data virt = dev_data.virt;

    float coef1 = dev_par.kernel.kernel_coefficient_1;
    float coef2 = dev_par.kernel.kernel_coefficient_2;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2f dx = x_j - x_i;

    float q_2 = dx.squaredNorm() * h_inv_2;
    float kernel_coef = (expf(-q_2) - coef1) * coef2 * coef3to2;

    Eigen::Vector2f vel_j = fluid.velocity[id_j];
    float pre_j = fluid.pressure[id_j];
    float rho_j = fluid.density[id_j];
    Eigen::Vector2f gravity(dev_par.physics.gravity[0], dev_par.physics.gravity[1]);

    sum_correction += kernel_coef;
    velocity += vel_j * kernel_coef;
    pressure += (pre_j - rho_j * gravity.dot(dx)) * kernel_coef;

    if (dx.squaredNorm() < nearest_distance_square)
    {
        virt.phase_type[id_i] = fluid.phase_type[id_j];
        nearest_distance_square = dx.squaredNorm();
    }
}

__global__ void getVirtData()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.virt.number)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Virtual_Data virt = dev_data.virt;
    Cell *cell_i = dev_data.cell + virt.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    Eigen::Vector2f velocity = Eigen::Vector2f::Zero();
    float pressure(0.0f), density(0.0f), sum_correction(0.0f);
    float nearest_distance_square(LARGE_FLOAT);

    Eigen::Vector2f x_i = virt.position[id];
    for (int n_id = 0; n_id < 9; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;
        int fluid_num_j = cell_j->fluid_num;
        for (int p_id = 0; p_id < fluid_num_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2f x_j = fluid.position[id_j];
            float dis_2 = (x_j - x_i).squaredNorm();
            if (dis_2 < impact_dis_2)
                virtDataCalculate(id, id_j, x_i, x_j, velocity, pressure, sum_correction, nearest_distance_square);
        }
    }

    if (sum_correction != 0.0f)
    {
        virt.velocity[id] = velocity / sum_correction;
        virt.pressure[id] = pressure / sum_correction;
    }
    else
    {
        virt.velocity[id] = Eigen::Vector2f::Zero();
        virt.pressure[id] = 0.0f;
    }

    switch (virt.boundary_type[id])
    {
    case SLIP:
        virt.velocity[id] = virt.velocity[id] - 2.0f * virt.velocity[id].dot(virt.normal[id]) * virt.normal[id];
        break;
    case NO_SLIP:
        virt.velocity[id] = -virt.velocity[id];
        break;
    }

    if (nearest_distance_square == 1.0e6)
        virt.phase_type[id] = LIQUID;
    virt.density[id] = pressureToDensity(virt.pressure[id], virt.phase_type[id]);
}

void DeviceCalculation::virtUpdate()
{
    if (!(Sim::parameter.hasVirtual()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int virt_num = Sim::parameter.virt.number;
    int virt_block = virt_num / thread_num + 1;
    getVirtData<<<virt_block, thread_num>>>();
    cudaDeviceSynchronize();
}