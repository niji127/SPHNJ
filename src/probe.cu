#include "probe.h"
#include "calculation.cuh"
#include "memory_copy.cuh"

#include <windows.h>

__device__ int findProbeCell(const Eigen::Vector2f &pos)
{
    float xmin = dev_par.domain.cell_xmin[0];
    float ymin = dev_par.domain.cell_xmin[1];
    float inter_inv_x = dev_par.domain.interval_inv[0];
    float inter_inv_y = dev_par.domain.interval_inv[1];
    int block_nx = dev_par.domain.cell_number[0];
    int block_ny = dev_par.domain.cell_number[1];

    return (int)((pos[1] - ymin) * inter_inv_y) + (int)((pos[0] - xmin) * inter_inv_x) * block_ny;
}

__global__ void getProbeData(Probe *probe)
{

    Eigen::Vector2f pos_i = probe->position;
    int cell_id = findProbeCell(pos_i);

    Cell *cell_i = dev_data.cell + cell_id;
    int fluid_num = cell_i->fluid_num;
    int solid_num = cell_i->solid_num;
    int virt_num = cell_i->virt_num;

    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Virtual_Data virt = dev_data.virt;

    float coef1 = dev_par.kernel.kernel_coefficient_1;
    float coef2 = dev_par.kernel.kernel_coefficient_2;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;
    float impact_dis_2 = dev_par.kernel.impact_distance_square;

    float sum_correction(0.0f);
    Eigen::Vector2f probe_velocity = Eigen::Vector2f::Zero();
    float probe_pressure(0.0f), probe_density(0.0f);
    float near_pos(10.0f);

    for (int n_id = 0; n_id < 9; n_id++)
    {
        int neighbor_cell_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_cell_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_cell_id;

        // fluid
        for (int i = 0; i < fluid_num; i++)
        {
            int id_j = cell_j->particle_list[i];
            Eigen::Vector2f pos_j = fluid.position[id_j];

            float dis_2 = (pos_j - pos_i).squaredNorm();
            if (dis_2 < impact_dis_2)
            {
                Eigen::Vector2f vel_j = fluid.velocity[id_j];
                float pre_j = fluid.pressure[id_j];
                float dens_j = fluid.density[id_j];
                float q_2 = dis_2 * h_inv_2;
                float kernel_coef = (expf(-q_2) - coef1) * coef2 * coef3to2;

                sum_correction += kernel_coef;
                probe_velocity += vel_j * kernel_coef;
                probe_pressure += pre_j * kernel_coef;
                probe_density += dens_j * kernel_coef;
            }
        }

        // solid
        for (int i = fluid_num; i < fluid_num + solid_num; i++)
        {
            int id_j = cell_j->particle_list[i];
            Eigen::Vector2f pos_j = solid.position[id_j].cast<float>();
            float dis_2 = (pos_j - pos_i).squaredNorm();
            if (dis_2 < impact_dis_2)
            {
                Eigen::Vector2f vel_j = solid.velocity[id_j].cast<float>();
                float pre_j = solid.pressure[id_j];
                float dens_j = 0.0f;
                float q_2 = dis_2 * h_inv_2;
                float kernel_coef = (expf(-q_2) - coef1) * coef2 * coef3to2;

                sum_correction += kernel_coef;
                probe_velocity += vel_j * kernel_coef;
                probe_pressure += pre_j * kernel_coef;
                probe_density += dens_j * kernel_coef;
            }
        }
    }
    if (sum_correction != 0.0f)
    {
        probe->velocity = probe_velocity / sum_correction;
        probe->pressure = probe_pressure / sum_correction;
        probe->density = probe_density / sum_correction;
    }
    else
    {
        probe->velocity = Eigen::Vector2f::Zero();
        probe->pressure = 0.0f;
        probe->density = dev_par.fluid.reference_density[0];
    }
}

void Probe::getData()
{
    Probe *device_probe(nullptr);
    cudaMalloc(&device_probe, sizeof(Probe));
    cudaMemcpy(device_probe, this, sizeof(Probe), cudaMemcpyHostToDevice);
    getProbeData<<<1, 1>>>(device_probe);
    cudaDeviceSynchronize();
    cudaMemcpy(this, device_probe, sizeof(Probe), cudaMemcpyDeviceToHost);
    cudaFree(device_probe);
}

void Probe::outputData()
{
    getData();
    std::string path = "..\\..\\out\\global";
    CreateDirectory(path.c_str(), NULL);
    std::fstream output;
    output.open(path + "\\probe-" + std::to_string(id) + ".dat", std::ios::app);
    output << position[0] << " " << position[1] << " " << velocity[0] << " " << velocity[1] << " " << pressure << " " << density << std::endl;
    output.close();
}