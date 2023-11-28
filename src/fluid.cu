#include "calculation.cuh"
#include "memory_copy.cuh"

#include "thrust/host_vector.h"

__device__ float pressureToDensity(float &pressure, const PhaseType &phase)
{
    float density_reference = dev_par.fluid.reference_density[phase];
    float coef = dev_par.fluid.coefficient_p2rho[phase];
    float gamma_inv = dev_par.fluid.gamma_inv[phase];

    return density_reference * powf(pressure * coef + 1.0f, gamma_inv);
}

__device__ float densityToPressure(float &density, const PhaseType &phase)
{
    float density_inverse = dev_par.fluid.reference_density_inverse[phase];
    float coef = dev_par.fluid.coefficient_rho2p[phase];
    float gamma = dev_par.fluid.gamma[phase];

    float pressure=(powf(density * density_inverse, gamma) - 1.0f) * coef;
    return max(pressure, dev_par.fluid.min_pressure);
}

__device__ void fluidFluidCalculate(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, float &velo_div, Eigen::Vector2f &du_dt, NearParticle &has_near_particle, float &color, Eigen::Vector2f &gradient, Eigen::Matrix2f &du_dx)
{
    Fluid_Data fluid = dev_data.fluid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2f u_i = fluid.velocity[id_i];
    Eigen::Vector2f u_j = fluid.velocity[id_j];

    float p_i = fluid.pressure[id_i];
    float p_j = fluid.pressure[id_j];
    PhaseType phase_i = fluid.phase_type[id_i];
    float cs = dev_par.fluid.sound_speed[phase_i];

    float rho_i = fluid.density[id_i];
    float rho_j = fluid.density[id_j];

    Eigen::Vector2f dx = x_j - x_i;
    float q_2 = dx.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = (u_i - u_j).dot(dx.normalized());
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * min(3.0f * du, cs);

    Eigen::Vector2f vs = (rho_i * u_i + rho_j * u_j) / (rho_i + rho_j);
    vs += (p_i - p_j) / (rho_i + rho_j) / cs * dx.normalized();

    du_dt += -2.0f * ps / rho_i * dx * kernel_coef;
    velo_div += 2.0f * rho_i * (u_i - vs).dot(dx) * kernel_coef;
    du_dx += (u_j - u_i) * dx.transpose() * kernel_coef;

    PhaseType phase_j = fluid.phase_type[id_j];
    if (!(phase_i == LIQUID && phase_j == GAS))
    {
        color += dx.dot(dx) * kernel_coef;
        gradient += dx * kernel_coef;
    }

    if (fluid.fluid_boundary_type[id_j] == FREE_SURFACE && phase_i == LIQUID && id_i != id_j)
    {
        if (dx.squaredNorm() < 0.25f * dev_par.kernel.impact_distance_square)
            has_near_particle = HAS_NEAR_PARTICLE;
        else if (has_near_particle == NO_PARTICLE)
            has_near_particle = HAS_PARTICLE;
    }

    if (!dev_par.fluid.is_fluid_viscid)
        return;
    float vis = (u_j - u_i).dot(dx);
    vis /= dx.squaredNorm() + 0.01f * dev_par.kernel.smoothing_length_square;
    du_dt += 8.0f * dev_par.fluid.viscosity[phase_i] * vis / rho_i * dx * kernel_coef;
}

__device__ void fluidSolidCalculate(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, float &velo_div, Eigen::Vector2f &du_dt, NearParticle &has_near_particle, float &color, Eigen::Vector2f &gradient, Eigen::Matrix2f &du_dx)
{
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2f u_i = fluid.velocity[id_i];
    Eigen::Vector2f u_j;

    Eigen::Vector2f normal = solid.normal[id_j];
    Eigen::Vector2f solid_vel = solid.velocity[id_j].cast<float>();
    Eigen::Vector2f normal_u = solid_vel.dot(normal) * normal;

    // u_j = u_i - 2.0f * (u_i - normal_u).dot(normal) * normal;
    u_j = solid_vel;

    Eigen::Vector2f dx = x_j - x_i;

    float p_j = solid.pressure[id_j];
    float p_i = fluid.pressure[id_i];
    PhaseType phase_i = fluid.phase_type[id_i];
    float cs = dev_par.fluid.sound_speed[phase_i];

    float rho_i = fluid.density[id_i];
    // float rho_j = pressureToDensity(p_j, phase_i);
    float rho_j = rho_i;

    float q_2 = dx.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = (u_i - u_j).dot(dx.normalized());
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * min(3.0f * du, cs);
    Eigen::Vector2f vs = (rho_i * u_i + rho_j * u_j) / (rho_i + rho_j);
    vs += (p_i - p_j) / (rho_i + rho_j) / cs * dx.normalized();

    du_dt += -2.0f * ps / rho_i * dx * kernel_coef;
    velo_div += 2.0f * rho_i * (u_i - vs).dot(dx) * kernel_coef;
    du_dx += (u_j - u_i) * dx.transpose() * kernel_coef;

    // color += dx.dot(dx) * kernel_coef;
    gradient += dx * kernel_coef;

    if (!dev_par.fluid.is_fluid_viscid)
        return;
    float vis = (u_j - u_i).dot(dx);
    vis /= dx.squaredNorm() + 0.01f * dev_par.kernel.smoothing_length_square;
    du_dt += 8.0f * dev_par.fluid.viscosity[phase_i] * vis / rho_i * dx * kernel_coef;
}

__device__ void fluidVirtCalculate(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, float &velo_div, Eigen::Vector2f &du_dt, NearParticle &has_near_particle, float &color, Eigen::Vector2f &gradient, Eigen::Matrix2f &du_dx)
{
    Fluid_Data fluid = dev_data.fluid;
    Virtual_Data virt = dev_data.virt;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;
    Eigen::Vector2f gravity = {dev_par.physics.gravity[0], dev_par.physics.gravity[1]};

    Eigen::Vector2f u_i = fluid.velocity[id_i];
    Eigen::Vector2f u_j = virt.velocity[id_j];

    Eigen::Vector2f dx = x_j - x_i;
    float q_2 = dx.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;

    float p_i = fluid.pressure[id_i];
    float p_j = virt.pressure[id_j];
    PhaseType phase_i = fluid.phase_type[id_i];
    float cs = dev_par.fluid.sound_speed[phase_i];

    float rho_i = fluid.density[id_i];
    float rho_j = dev_par.fluid.reference_density[virt.phase_type[id_j]];

    // if (x_j(1) > 2.0f)
    // {
    //     u_j = u_i;
    //     p_j = 0.0;
    // }

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = (u_i - u_j).dot(dx.normalized());
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * min(3.0f * du, cs);

    Eigen::Vector2f vs = (rho_i * u_i + rho_j * u_j) / (rho_i + rho_j);
    vs += (p_i - p_j) / (rho_i + rho_j) / cs * dx.normalized();

    du_dt += -2.0f * ps / rho_i * dx * kernel_coef;
    velo_div += 2.0f * rho_i * (u_i - vs).dot(dx) * kernel_coef;
    du_dx += (u_j - u_i) * dx.transpose() * kernel_coef;

    PhaseType phase_j = fluid.phase_type[id_j];
    if (!(phase_i == LIQUID && phase_j == GAS))
    {
        color += dx.dot(dx) * kernel_coef;
        gradient += dx * kernel_coef;
    }

    if (!dev_par.fluid.is_fluid_viscid)
        return;
    float vis = (u_j - u_i).dot(dx);
    vis /= dx.squaredNorm() + 0.01f * dev_par.kernel.smoothing_length_square;
    du_dt += 8.0f * dev_par.fluid.viscosity[phase_i] * vis / rho_i * dx * kernel_coef;
}

__global__ void findFluidPair()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.fluid.number_total)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Virtual_Data virt = dev_data.virt;
    Cell *cell_i = dev_data.cell + fluid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    float velo_div(0.0f), color(0.0f);
    Eigen::Vector2f gradient = Eigen::Vector2f::Zero();
    Eigen::Vector2f du_dt = Eigen::Vector2f::Zero();
    Eigen::Matrix2f du_dx = Eigen::Matrix2f::Zero();
    NearParticle has_near_particle = NO_PARTICLE;
    Eigen::Vector2f x_i = fluid.position[id];
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
                fluidFluidCalculate(id, id_j, x_i, x_j, velo_div, du_dt, has_near_particle, color, gradient, du_dx);
        }

        int solid_num_j = cell_j->solid_num;
        for (int p_id = fluid_num_j; p_id < fluid_num_j + solid_num_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2f x_j = solid.position[id_j].cast<float>();
            float dis_2 = (x_j - x_i).squaredNorm();
            if (dis_2 < impact_dis_2)
                fluidSolidCalculate(id, id_j, x_i, x_j, velo_div, du_dt, has_near_particle, color, gradient, du_dx);
        }

        int virt_sum_j = cell_j->virt_num;
        for (int p_id = fluid_num_j + solid_num_j; p_id < fluid_num_j + solid_num_j + virt_sum_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2f x_j = virt.position[id_j];
            float dis_2 = (x_j - x_i).squaredNorm();
            if (dis_2 < impact_dis_2)
                fluidVirtCalculate(id, id_j, x_i, x_j, velo_div, du_dt, has_near_particle, color, gradient, du_dx);
        }
    }

    fluid.drho_dt[id] = velo_div;
    fluid.du_dt[id] = du_dt;
    fluid.has_near_particle[id] = has_near_particle;
    fluid.color[id] = color;
    fluid.gradient[id] = gradient;
    fluid.du_dx[id] = du_dx;
}

__global__ void fluidUpdate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.fluid.number_total)
        return;

    Fluid_Data fluid = dev_data.fluid;
    PhaseType phase = fluid.phase_type[id];
    Eigen::Vector2f gravity = Eigen::Vector2f(dev_par.physics.gravity[0], dev_par.physics.gravity[1]);
    Eigen::Vector2f xmin = Eigen::Vector2f(dev_par.domain.domain_xmin[0], dev_par.domain.domain_xmin[1]);
    Eigen::Vector2f xmax = Eigen::Vector2f(dev_par.domain.domain_xmax[0], dev_par.domain.domain_xmax[1]);

    float dt = dev_par.time.dt;

    if (dev_par.physics.trasient_type == STEADY)
        fluid.du_dt[id] += -0.1f * fluid.velocity[id];
    fluid.position[id] += fluid.velocity[id] * dt;
    fluid.velocity[id] += (fluid.du_dt[id] + gravity) * dt;
    fluid.density[id] += fluid.drho_dt[id] * dt;
    fluid.density[id] = max(fluid.density[id], dev_par.fluid.min_density[phase]);
    fluid.pressure[id] = densityToPressure(fluid.density[id], phase);

    // saturate
    fluid.position[id](0) = max(fluid.position[id](0), xmin[0]);
    fluid.position[id](0) = min(fluid.position[id](0), xmax[0]);
    fluid.position[id](1) = max(fluid.position[id](1), xmin[1]);
    fluid.position[id](1) = min(fluid.position[id](1), xmax[1]);

    // float max_mach_number = 2.0;
    // if (fluid.velocity[id].norm() > max_mach_number * dev_par.fluid.sound_speed[phase])
    //     fluid.velocity[id] *= dev_par.fluid.sound_speed[phase] / fluid.velocity[id].norm();
}

void DeviceCalculation::fluidCalculation()
{
    if (!(Sim::parameter.hasFluid()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int fluid_num = Sim::parameter.fluid.number_total;
    int fluid_block = fluid_num / thread_num + 1;

    findFluidPair<<<fluid_block, thread_num>>>();
    cudaDeviceSynchronize();
    fluidUpdate<<<fluid_block, thread_num>>>();
    cudaDeviceSynchronize();
}

__device__ void fluidCorrectionCal(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, Eigen::Vector2f &gradient)
{
    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2f dx = x_j - x_i;
    float q_2 = dx.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;
    gradient += dx * kernel_coef;
}

__global__ void correctionCalculate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.fluid.number_total)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Virtual_Data virt = dev_data.virt;
    Cell *cell_i = dev_data.cell + fluid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    Eigen::Vector2f gradient = Eigen::Vector2f::Zero();
    Eigen::Vector2f x_i = fluid.position[id];
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
                fluidCorrectionCal(id, id_j, x_i, x_j, gradient);
        }

        int solid_num_j = cell_j->solid_num;
        for (int p_id = fluid_num_j; p_id < fluid_num_j + solid_num_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2f x_j = solid.position[id_j].cast<float>();
            float dis_2 = (x_j - x_i).squaredNorm();
            if (dis_2 < impact_dis_2)
                fluidCorrectionCal(id, id_j, x_i, x_j, gradient);
        }

        int virt_sum_j = cell_j->virt_num;
        for (int p_id = fluid_num_j + solid_num_j; p_id < fluid_num_j + solid_num_j + virt_sum_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2f x_j = virt.position[id_j];
            float dis_2 = (x_j - x_i).squaredNorm();
            if (dis_2 < impact_dis_2)
                fluidCorrectionCal(id, id_j, x_i, x_j, gradient);
        }
    }
    fluid.gradient[id] = gradient;
}

__device__ Eigen::Vector2f shiftDisplacementCalculate(int id)
{
    Fluid_Data fluid = dev_data.fluid;

    float coef = dev_par.fluid.particle_shift_coef;
    float hsml = dev_par.kernel.smoothing_length;
    float dt = dev_par.time.dt;
    Eigen::Vector2f displacement = -coef * hsml * dt * fluid.gradient[id];

    if (fluid.phase_type[id] == GAS)
    {
        fluid.fluid_boundary_type[id] = INNER;
        return displacement;
    }

    if (fluid.color[id] < 1.5f)
    {
        if (fluid.has_near_particle[id] == NO_PARTICLE)
        {
            fluid.fluid_boundary_type[id] = SPLASH;
            return Eigen::Vector2f::Zero();
        }
        else
        {
            fluid.fluid_boundary_type[id] = FREE_SURFACE;
            Eigen::Vector2f normal = fluid.gradient[id].normalized();
            if (normal.norm() < EPS_FOR_SUM)
                normal.setZero();
            return (Eigen::Matrix2f::Identity() - normal * normal.transpose()) * displacement;
        }
    }

    if (fluid.color[id] < 2.0f)
    {
        if (fluid.has_near_particle[id] == HAS_NEAR_PARTICLE)
        {
            fluid.fluid_boundary_type[id] = VICINITY;
            Eigen::Vector2f normal = fluid.gradient[id].normalized();
            if (normal.norm() < EPS_FOR_SUM)
                normal.setZero();
            return (Eigen::Matrix2f::Identity() - normal * normal.transpose()) * displacement;
        }
    }
    fluid.fluid_boundary_type[id] = INNER;
    return displacement;
}

__global__ void correctionUpdate(float max_speed)
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.fluid.number_total)
        return;

    Fluid_Data fluid = dev_data.fluid;

    Eigen::Vector2f shift_displacement = max_speed * shiftDisplacementCalculate(id);
    fluid.position[id] += shift_displacement;
    fluid.velocity[id] += fluid.du_dx[id] * shift_displacement;
}

float getMaxSpeed(Sim *sim)
{
    thrust::host_vector<Eigen::Vector2f> host_velocity = sim->fluid_list.velocity;
    float max_speed = 0.0f;
    for (auto fp : host_velocity)
    {
        if (fp.norm() > max_speed)
            max_speed = fp.norm();
    }
    float min_sound_speed = min(Sim::parameter.fluid.sound_speed[0], Sim::parameter.fluid.sound_speed[1]);
    return min(max_speed, 0.1f * min_sound_speed);
}

void DeviceCalculation::fluidCorrection(Sim *sim)
{
    if (!(Sim::parameter.hasFluid()))
        return;

    float max_speed = getMaxSpeed(sim);

    int thread_num = Sim::parameter.kernel.thread_num;
    int fluid_num = Sim::parameter.fluid.number_total;
    int fluid_block = fluid_num / thread_num + 1;
    for (int i = 0; i < Sim::parameter.time.shift_per_step; i++)
    {
        correctionCalculate<<<fluid_block, thread_num>>>();
        cudaDeviceSynchronize();
        correctionUpdate<<<fluid_block, thread_num>>>(max_speed);
        cudaDeviceSynchronize();
    }
}