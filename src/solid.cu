#include "calculation.cuh"
#include "memory_copy.cuh"

extern __device__ float densityToPressure(float &density, const PhaseType &phase);
extern __device__ float pressureToDensity(float &pressure, const PhaseType &phase);

__device__ void solidSolidDeformation(int id_i, int id_j, const Eigen::Vector2d &x_ref_i, const Eigen::Vector2d &x_ref_j, Eigen::Matrix2d &deformation, Eigen::Matrix2d &correction)
{
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2d dx = solid.position[id_j] - solid.position[id_i];

    Eigen::Vector2d dx_ref = x_ref_j - x_ref_i;
    float q_2 = dx_ref.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;

    deformation += dx * dx_ref.transpose() * kernel_coef;
    correction += dx_ref * dx_ref.transpose() * kernel_coef;
}

__device__ float getVonMisesStressSquare(const Eigen::Matrix2d &stress)
{
    float mises = 0.0f;
    mises += 0.5f * powf(stress(0, 0) - stress(1, 1), 2.0f);
    mises += 3.0f * powf(stress(0, 1), 2.0f);
    return mises;
}

__device__ float getVonMisesStress(const Eigen::Matrix2d &stress)
{
    float mises = 0.0f;
    mises += 0.5f * powf(stress(0, 0) - stress(1, 1), 2.0f);
    mises += 3.0f * powf(stress(0, 1), 2.0f);
    return sqrtf(mises);
}

__device__ double hyperElasticFunction(double strain, const Eigen::Matrix2d &deformation)
{
    float yield_stress = dev_par.solid.yield_stress;

    double true_strain = sqrt(1.0 + 2.0 * strain) - 1.0;
    if (true_strain > 0.0f)
    {
        float yield_strain = yield_stress / dev_par.solid.auxetic_modulus;
        if (true_strain < yield_strain)
        {
            return true_strain * dev_par.solid.auxetic_modulus;
        }
        float harden_strain = dev_par.solid.harden_strain;
        float yield_modulus = dev_par.solid.auxetic_modulus * 0.0f;
        if (true_strain < harden_strain)
        {
            return dev_par.solid.yield_stress + yield_modulus * (true_strain - yield_strain);
        }
        float harden_coef = dev_par.solid.harden_coef;
        return harden_coef * yield_stress * pow(true_strain - harden_strain, 2.0) + yield_stress + yield_modulus * (harden_strain - yield_strain);
    }
    else
    {
        float yield_strain = -yield_stress / dev_par.solid.auxetic_modulus;
        if (true_strain > yield_strain)
        {
            return true_strain * dev_par.solid.auxetic_modulus;
        }
        float harden_strain = -dev_par.solid.harden_strain;
        float yield_modulus = dev_par.solid.auxetic_modulus * 0.0f;
        if (true_strain > harden_strain)
        {
            return -(dev_par.solid.yield_stress + yield_modulus * (yield_strain - true_strain));
        }
        float harden_coef = dev_par.solid.harden_coef;
        return -(harden_coef * yield_stress * pow(true_strain - harden_strain, 2.0) + yield_stress + yield_modulus * (yield_strain - harden_strain));
    }
}

__device__ double plasticElasticFunction(double strain, const Eigen::Matrix2d &deformation, double &max_strain)
{
    float modulus = dev_par.solid.auxetic_modulus;
    float yield_stress = dev_par.solid.yield_stress;
    float yield_strain = dev_par.solid.yield_strain;
    float harden_strain = dev_par.solid.harden_strain;
    float harden_coef = dev_par.solid.harden_coef;
    double true_strain;
    if (strain > -0.5)
        true_strain = sqrt(1.0 + 2.0 * strain) - 1.0;
    else
        true_strain = -1.0;
    if (true_strain > max_strain)
    {
        return modulus * (true_strain - max_strain) - yield_stress;
    }
    else if (true_strain > -harden_strain)
    {
        if (max_strain > true_strain && dev_par.solid.plastic_type == PLASTIC)
            max_strain = true_strain;
        return -yield_stress;
    }
    else
    {
        if (max_strain > true_strain && dev_par.solid.plastic_type == PLASTIC)
            max_strain = true_strain;
        return -(harden_coef * pow(fabs(-harden_strain - true_strain), 2.0) + yield_stress);
    }
}

__device__ Eigen::Matrix2d getMatrixEigenValues(const Eigen::Matrix2d &matrix)
{
    Eigen::Matrix2d eigen_value = Eigen::Matrix2d::Zero();
    eigen_value(0, 0) = 0.5 * (matrix(0, 0) + matrix(1, 1));
    eigen_value(0, 0) += sqrt(0.25 * pow(matrix(0, 0) - matrix(1, 1), 2.0) + matrix(0, 1) * matrix(1, 0));
    eigen_value(1, 1) = 0.5 * (matrix(0, 0) + matrix(1, 1));
    eigen_value(1, 1) -= sqrt(0.25 * pow(matrix(0, 0) - matrix(1, 1), 2.0) + matrix(0, 1) * matrix(1, 0));
    return eigen_value;
}

__device__ Eigen::Matrix2d getMatrixEigenVectors(const Eigen::Matrix2d &matrix)
{
    double alpha = 0.5 * atan2(matrix(0, 1) + matrix(1, 0), matrix(0, 0) - matrix(1, 1));
    Eigen::Matrix2d vectors;
    vectors(0, 0) = cos(alpha);
    vectors(0, 1) = -sin(alpha);
    vectors(1, 0) = sin(alpha);
    vectors(1, 1) = cos(alpha);
    return vectors;
}

__device__ Eigen::Matrix2d auxeticStressCalculation(const Eigen::Matrix2d &deformation, Eigen::Matrix2d &max_strain)
{
    Eigen::Matrix2d green_strain = 0.5 * (deformation.transpose() * deformation - Eigen::Matrix2d::Identity());

    Eigen::Matrix2d principle_strain = getMatrixEigenValues(green_strain);
    Eigen::Matrix2d rotation_matrix = getMatrixEigenVectors(green_strain);
    Eigen::Matrix2d max_strain_principle = rotation_matrix.inverse() * max_strain * rotation_matrix;

    Eigen::Matrix2d principle_stress = Eigen::Matrix2d::Zero();
    for (int i = 0; i < 2; i++)
    {
        // principle_stress(i, i) = hyperElasticFunction(principle_strain(i, i), deformation);
        principle_stress(i, i) = plasticElasticFunction(principle_strain(i, i), deformation, max_strain_principle(i, i));
    }
    max_strain = rotation_matrix * max_strain_principle * rotation_matrix.inverse();
    // principle_stress /= deformation.determinant();
    return deformation * rotation_matrix * principle_stress * rotation_matrix.inverse();
}

__device__ Eigen::Matrix2d linearStressCalculation(const Eigen::Matrix2d &deformation)
{
    float lambda = dev_par.solid.lambda;
    float nu = dev_par.solid.nu;

    Eigen::Matrix2d green_strain = 0.5 * (deformation.transpose() * deformation - Eigen::Matrix2d::Identity());
    Eigen::Matrix2d stress_2nd_piola = 2.0 * nu * green_strain + lambda * green_strain.trace() * Eigen::Matrix2d::Identity();
    return deformation * stress_2nd_piola;
}

__device__ Eigen::Matrix2d solidStressCalculation(const Eigen::Matrix2d &deformation, Eigen::Matrix2d &max_strain, const MaterialType &material_type)
{
    if (material_type == AUXETIC)
        return auxeticStressCalculation(deformation, max_strain);
    else if (material_type == LINEAR_ELASTIC)
        return linearStressCalculation(deformation);
    return Eigen::Matrix2d::Zero();
}

__global__ void solidDeformationCalculate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;

    Solid_Data solid = dev_data.solid;
    Cell *cell_i = dev_data.cell + solid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    Eigen::Matrix2d correction = Eigen::Matrix2d::Zero();
    Eigen::Matrix2d deformation = Eigen::Matrix2d::Zero();

    Eigen::Vector2d x_ref_i = solid.reference_position[id];
    // deformation && correction
    for (int n_id = 0; n_id < 9; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;
        int fluid_num_j = cell_j->fluid_num;
        int solid_num_j = cell_j->solid_num;
        for (int p_id = fluid_num_j; p_id < fluid_num_j + solid_num_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2d x_ref_j = solid.reference_position[id_j];
            float dis_ref_2 = (x_ref_j - x_ref_i).squaredNorm();
            if (dis_ref_2 < impact_dis_2)
                solidSolidDeformation(id, id_j, x_ref_i, x_ref_j, deformation, correction);
        }
    }

    if (correction.determinant() != 0.0f)
        correction = correction.inverse();
    else
        correction = Eigen::Matrix2d::Identity();
    solid.deformation[id] = deformation * correction;
    solid.stress_1st_piola[id] = solidStressCalculation(deformation * correction, solid.max_strain[id], solid.material_type[id]);
    solid.correction[id] = correction;
}

__device__ void solidFluidCouple(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, Eigen::Vector2f &du_dt)
{
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2f u_i;
    Eigen::Vector2f u_j = fluid.velocity[id_j];

    Eigen::Vector2f normal = solid.normal[id_i];
    Eigen::Vector2f solid_vel = solid.velocity[id_i].cast<float>();
    Eigen::Vector2f normal_u = solid_vel.dot(normal) * normal;

    // u_i = u_j - 2.0f * (u_j - normal_u).dot(normal) * normal;
    u_i = solid_vel;

    Eigen::Vector2f dx = x_j - x_i;
    float q_2 = dx.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;

    float p_i = solid.pressure[id_i];
    float p_j = fluid.pressure[id_j];
    PhaseType phase = fluid.phase_type[id_j];
    float cs = dev_par.fluid.sound_speed[phase];

    // float rho_i = pressureToDensity(p_i, phase);
    float rho_j = fluid.density[id_j];
    float rho_i = rho_j;

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = (u_i - u_j).dot(dx.normalized());
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * min(3.0f * du, cs);

    float rho_ref;
    if (solid.material_type[id_i] == LINEAR_ELASTIC)
        rho_ref = dev_par.solid.reference_density;
    else if (solid.material_type[id_i] == AUXETIC)
        rho_ref = dev_par.solid.auxetic_density;

    du_dt += -2.0f * ps / rho_ref * dx * kernel_coef;
}

__device__ void solidSolidForce(int id_i, int id_j, Eigen::Vector2d &x_ref_i, Eigen::Vector2d &x_ref_j, Eigen::Vector2d &du_dt)
{
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;
    float rho_ref;
    if (solid.material_type[id_i] == LINEAR_ELASTIC)
        rho_ref = dev_par.solid.reference_density;
    else if (solid.material_type[id_i] == AUXETIC)
        rho_ref = dev_par.solid.auxetic_density;

    Eigen::Vector2d dx_ref = x_ref_j - x_ref_i;
    float q_2 = dx_ref.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;

    Eigen::Matrix2d deformation_i = solid.deformation[id_i];
    Eigen::Matrix2d stress_i = solid.stress_1st_piola[id_i];
    Eigen::Matrix2d stress_j = solid.stress_1st_piola[id_j];
    Eigen::Matrix2d correction_i = solid.correction[id_i];
    Eigen::Matrix2d correction_j = solid.correction[id_j];

    du_dt += (stress_j * correction_j + stress_i * correction_i) * dx_ref * kernel_coef / rho_ref;
    // if (dx_ref(1) > 0.0)
    // {
    //     solid.couple_dudt[id_i] += ((stress_j * correction_j + stress_i * correction_i) * dx_ref * kernel_coef / rho_ref).cast<float>();
    // }

    Eigen::Vector2d dx = solid.position[id_j] - solid.position[id_i];
    Eigen::Vector2d du = solid.velocity[id_j] - solid.velocity[id_i];

    float du_dot_dx = du.dot(dx);
    // if (du_dot_dx > 0.0f)
    //     return;

    float alpha = dev_par.solid.artificial_viscocity[0];
    float beta = dev_par.solid.artificial_viscocity[1];
    float eps = 0.01f * dev_par.kernel.smoothing_length_square;
    float cs = dev_par.solid.sound_speed;

    float dis_square = dx.squaredNorm();
    float nu = dev_par.kernel.smoothing_length * du_dot_dx / (dis_square + eps);

    float viscosity = alpha * nu * cs;
    viscosity -= beta * nu * nu;

    du_dt += deformation_i.adjoint().transpose() * viscosity * correction_i * dx_ref * kernel_coef;
}

__global__ void solidForceCalculate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Cell *cell_i = dev_data.cell + solid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;

    Eigen::Vector2f couple_du_dt = Eigen::Vector2f::Zero();
    Eigen::Vector2d du_dt = Eigen::Vector2d::Zero();

    Eigen::Vector2f x_i = solid.position[id].cast<float>();
    Eigen::Vector2d x_ref_i = solid.reference_position[id];

    // du_dt
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
                solidFluidCouple(id, id_j, x_i, x_j, couple_du_dt);
        }
        int solid_num_j = cell_j->solid_num;
        for (int p_id = fluid_num_j; p_id < fluid_num_j + solid_num_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2d x_ref_j = solid.reference_position[id_j];
            float dis_ref_2 = (x_ref_j - x_ref_i).squaredNorm();
            if (dis_ref_2 < impact_dis_2)
                solidSolidForce(id, id_j, x_ref_i, x_ref_j, du_dt);
        }
    }
    solid.couple_dudt[id] = couple_du_dt;
    solid.du_dt[id] = du_dt + couple_du_dt.cast<double>();
}

__device__ void solidNormalCalculate(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, Eigen::Vector2f &normal)
{
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2f dx = x_j - x_i;

    float q_2 = dx.squaredNorm() * h_inv_2;
    float kernel_coef = expf(-q_2) * coef * coef3to2;

    normal -= dx * kernel_coef;
}

__device__ void solidPressureCalculate(int id_i, int id_j, Eigen::Vector2f &x_i, Eigen::Vector2f &x_j, float &pressure, float &sum_correction)
{
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;

    float coef1 = dev_par.kernel.kernel_coefficient_1;
    float coef2 = dev_par.kernel.kernel_coefficient_2;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float coef3to2 = dev_par.kernel.coefficient_3d_to_2d;

    Eigen::Vector2f dx = x_j - x_i;

    float q_2 = dx.squaredNorm() * h_inv_2;
    float pre_j = fluid.pressure[id_j];
    float rho_j = fluid.density[id_j];
    float kernel_coef = (expf(-q_2) - coef1) * coef2 * coef3to2;

    Eigen::Vector2f accel = solid.du_dt[id_i].cast<float>();

    sum_correction += kernel_coef;
    // pressure += (pre_j + rho_j * accel.dot(dx)) * kernel_coef;
    pressure += pre_j * kernel_coef;
}

__global__ void solidCoupleUpdate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Cell *cell_i = dev_data.cell + solid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    Eigen::Vector2f normal = Eigen::Vector2f::Zero();
    float pressure(0.0f), sum_correction(0.0f);

    Eigen::Vector2f x_i = solid.position[id].cast<float>();
    // normal
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
                solidPressureCalculate(id, id_j, x_i, x_j, pressure, sum_correction);
        }
        int solid_num_j = cell_j->solid_num;
        for (int p_id = fluid_num_j; p_id < fluid_num_j + solid_num_j; p_id++)
        {
            int id_j = cell_j->particle_list[p_id];
            Eigen::Vector2f x_j = solid.position[id_j].cast<float>();
            float dis_2 = (x_j - x_i).squaredNorm();
            if (dis_2 < impact_dis_2)
                solidNormalCalculate(id, id_j, x_i, x_j, normal);
        }
    }

    if (normal.norm() < EPS_FOR_SUM)
        solid.normal[id] = Eigen::Vector2f::Zero();
    else
        solid.normal[id] = normal.normalized();
    if (sum_correction != 0.0f)
        solid.pressure[id] = max(pressure / sum_correction, dev_par.fluid.min_pressure);
    else
        solid.pressure[id] = 0.0f;
}

__global__ void solidUpdate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;
    Solid_Data solid = dev_data.solid;
    Eigen::Vector2f gravity = Eigen::Vector2f(dev_par.physics.gravity[0], dev_par.physics.gravity[1]);
    Eigen::Vector2f xmin = Eigen::Vector2f(dev_par.domain.domain_xmin[0], dev_par.domain.domain_xmin[1]);
    Eigen::Vector2f xmax = Eigen::Vector2f(dev_par.domain.domain_xmax[0], dev_par.domain.domain_xmax[1]);

    float dt = dev_par.time.dt / (float)dev_par.time.solid_sub_step;
    switch (solid.boundary_type[id])
    {
    case FIXED_PARTICLE:
        solid.velocity[id] = Eigen::Vector2d(dev_par.solid.fixed_speed[0], dev_par.solid.fixed_speed[1]);
        solid.position[id] += solid.velocity[id] * dt;
        solid.du_dt[id] = Eigen::Vector2d::Zero();
        break;
    case UNCONSTRAINED_PARTICLE:
        if (dev_par.physics.trasient_type == STEADY)
            solid.du_dt[id] += -0.1 * solid.velocity[id];
        solid.position[id] += solid.velocity[id] * dt;
        solid.velocity[id] += (solid.du_dt[id] + gravity.cast<double>()) * dt;
        break;
    }
}

void DeviceCalculation::solidCalculation(Sim *sim)
{
    if (!(Sim::parameter.hasSolid()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int solid_num = Sim::parameter.solid.number;
    int solid_block = solid_num / thread_num + 1;

    for (int i = 0; i < Sim::parameter.time.solid_sub_step; i++)
    {
        SimTime *time = Sim::parameter.getTime();
        solidDeformationCalculate<<<solid_block, thread_num>>>();
        cudaDeviceSynchronize();
        solidForceCalculate<<<solid_block, thread_num>>>();
        cudaDeviceSynchronize();
        solidUpdate<<<solid_block, thread_num>>>();
        cudaDeviceSynchronize();
    }
}

void DeviceCalculation::coupleCalculation()
{
    if (!(Sim::parameter.hasFluid()))
        return;
    if (!(Sim::parameter.hasSolid()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int solid_num = Sim::parameter.solid.number;
    int solid_block = solid_num / thread_num + 1;

    solidCoupleUpdate<<<solid_block, thread_num>>>();
    cudaDeviceSynchronize();
}