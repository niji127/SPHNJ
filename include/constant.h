#include "parameter.h"
#include "particle.h"

struct PhysicsConstant
{
    float gravity[2];
};

struct FluidConstant
{
    int number[2];
    int number_total;
    float sound_speed[2], reference_density[2];
    float reference_density_inverse[2];
    float gamma[2], viscosity[2]; // dynamic viscosity
    float particle_shift_coef;
    float gamma_inv[2];
    float coefficient_p2rho[2];
    float coefficient_rho2p[2];
    float min_pressure_for_gas;
    bool is_fluid_viscid;
};

struct SolidConstant
{
    int number;
    float density, modulus, poisson_ratio;
    float artificial_viscocity[2];
    bool is_nonlinear;
};

struct VirtConstant
{
    int number;
};

struct Domain
{
    float xmin[2], xmax[2];
    float boundary_xmin[2], boundary_xmax[2];
    int block_number[2];
    float interval[2];
    float interval_inv[2];
    int block_number_total;
};

struct TimeConstant
{
    float dt, current_time;
    int istart, iend, i;
    int file_step;
    int solid_sub_step, shift_sub_step;
    int result_per_step;
}

class Kernel
{
public:
    int solid_fluid_radius_ratio;
    float fluid_smoothing_length;
    float fluid_particle_diameter;
    float fluid_particle_volume;
    float fluid_diameter_square;

    float solid_smoothing_length;
    float solid_particle_diameter;
    float solid_particle_volume;
    float solid_diameter_square;

    int block_size;
    int fluid_pair_size;
    int solid_pair_size;
    int thread_num;

    float fluid_smoothing_length_inverse;
    float solid_smoothing_length_inverse;

    float fluid_smoothing_length_inverse_square;
    float solid_smoothing_length_inverse_square;

    float fluid_impact_distance_square;
    float solid_impact_distance_square;

    float kernel_coefficient_1;
    float kernel_coefficient_2;

    float kernel_differential_fluid;
    float kernel_differential_solid;
    float coefficient_3d_to_2d;

    bool readParameter(const Config &config);
    void initiate();
    bool readNumber(const std::string &path);
};

    class CudaConstant
    {
    public:
        int solid_fluid_radius_ratio;
        float fluid_smoothing_length;
        float fluid_particle_diameter;
        float fluid_particle_volume;
        float fluid_diameter_square;

        float solid_smoothing_length;
        float solid_particle_diameter;
        float solid_particle_volume;
        float solid_diameter_square;

        int block_size;
        int fluid_pair_size;
        int solid_pair_size;
        int thread_num;

        float fluid_smoothing_length_inverse;
        float solid_smoothing_length_inverse;
        float fluid_smoothing_length_inverse_square;
        float solid_smoothing_length_inverse_square;

        float fluid_impact_distance_square;
        float solid_impact_distance_square;

        float kernel_coefficient_1;
        float kernel_coefficient_2;

        float kernel_differential_fluid;
        float kernel_differential_solid;
        float coefficient_3d_to_2d;

        Fluid_Particle *fluid_particle;
        Solid_Particle *solid_particle;
        Virtual_Particle *virtual_particle;

        float xmin, xmax, ymin, ymax;
        float gravity[2];

        void copyFromPhysics(Physics *physics)
        {
            for (int i = 0; i < 2; i++)
                this->gravity[i] = physics->gravity[0];
        }

        void copyFromDomain(Domain *domain)
        {
            this->xmin = domain->boundary_xmin[0];
            this->xmax = domain->boundary_xmax[0];
            this->ymin = domain->boundary_xmin[1];
            this->ymax = domain->boundary_xmax[1];
        }

        void copyFromKernel(Kernel *kernel)
        {
            this->solid_fluid_radius_ratio = kernel->solid_fluid_radius_ratio;
            this->fluid_smoothing_length = kernel->fluid_smoothing_length;
            this->fluid_particle_diameter = kernel->fluid_particle_diameter;
            this->fluid_particle_volume = kernel->fluid_particle_volume;
            this->fluid_diameter_square = kernel->fluid_diameter_square;

            this->solid_smoothing_length = kernel->solid_smoothing_length;
            this->solid_particle_diameter = kernel->solid_particle_diameter;
            this->solid_particle_volume = kernel->solid_particle_volume;
            this->solid_diameter_square = kernel->solid_diameter_square;

            this->block_size = kernel->block_size;
            this->fluid_pair_size = kernel->fluid_pair_size;
            this->solid_pair_size = kernel->solid_pair_size;
            this->thread_num = kernel->thread_num;

            this->fluid_smoothing_length_inverse = kernel->fluid_smoothing_length_inverse;
            this->solid_smoothing_length_inverse = kernel->solid_smoothing_length_inverse;

            this->fluid_smoothing_length_inverse_square = kernel->fluid_smoothing_length_inverse_square;
            this->solid_smoothing_length_inverse_square = kernel->solid_smoothing_length_inverse_square;

            this->fluid_impact_distance_square = kernel->fluid_impact_distance_square;
            this->solid_impact_distance_square = kernel->solid_impact_distance_square;

            this->kernel_coefficient_1 = kernel->kernel_coefficient_1;
            this->kernel_coefficient_2 = kernel->kernel_coefficient_2;

            this->kernel_differential_fluid = kernel->kernel_differential_fluid;
            this->kernel_differential_solid = kernel->kernel_differential_solid;
            this->coefficient_3d_to_2d = kernel->coefficient_3d_to_2d;
        }
    };