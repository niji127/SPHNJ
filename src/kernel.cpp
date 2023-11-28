#include "kernel.h"
#include "global.h"

bool Kernel::readNumber(const std::string &path)
{
    return false;
}

bool Kernel::readParameter(const Config &config)
{
    smoothing_length = config.Read("smoothing_length", 0.0f);
    particle_diameter = config.Read("particle_diameter", 0.0f);
    impact_length_hsml_ratio = config.Read("impact_length_hsml_ratio", 3.0f);
    thread_num = config.Read("thread_num", 32);

    if (smoothing_length == 0.0f || particle_diameter == 0.0f)
    {
        std::cout << "failed to read kernel information" << std::endl;
        return true;
    }

    return false;
}

void Kernel::initiate()
{
    particle_volume = powf(particle_diameter, 3.0f);
    diameter_square = powf(particle_diameter, 2.0f);

    smoothing_length_square = powf(smoothing_length, 2.0f);
    smoothing_length_inverse = 1.0f / smoothing_length;
    smoothing_length_inverse_square = powf(smoothing_length_inverse, 2.0f);
    impact_distance_square = powf(impact_length_hsml_ratio, 2.0f) * pow(smoothing_length, 2.0f);

    float hsml_diameter_ratio = smoothing_length / particle_diameter;

    kernel_coefficient_1 = expf(-powf(impact_length_hsml_ratio, 2.0f));
    kernel_coefficient_2 = 1.0f / tg_coef_2 * pow(hsml_diameter_ratio, -3.0f);
    kernel_differential = 2.0f / tg_coef_2 * pow(hsml_diameter_ratio, -3.0f) * pow(smoothing_length, -2.0f);
    coefficient_3d_to_2d = COEF_3D_TO_2D;
}
