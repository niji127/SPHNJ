#include "properties.h"
#include <iostream>
#include <fstream>
#include <string>
#include <windows.h>

bool FluidProperties::readNumber(const std::string &path)
{
    std::ifstream infile;
    infile.open(path + "\\fluid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read fluid particles parameter" << std::endl;
        return true;
    }
    std::string str;
    std::getline(infile, str);
    std::istringstream string_to_num(str);
    string_to_num >> number[0] >> number[1];
    number_total = number[0] + number[1];
    if (number_total == 0)
        return false;
    int lines_check(0);
    while (std::getline(infile, str))
        lines_check++;
    if (lines_check != number_total)
    {
        std::cout << "mismatch fluid particle numbers" << std::endl;
        return true;
    }
    infile.close();
    return false;
}

bool FluidProperties::readParameter(const Config &config)
{
    particle_shift_coef = config.Read("particle_shift_coef", 0.4f);
    is_fluid_viscid = config.Read("is_viscid", false);
    min_pressure = config.Read("min_pressure", -1.013e5f);

    sound_speed[0] = config.Read("liquid_sound_speed", 1500.0f);
    reference_density[0] = config.Read("liquid_density", 1000.0f);
    gamma[0] = config.Read("liquid_gamma", 7.0f);
    viscosity[0] = config.Read("liquid_viscosity", 1.01e-3f); // 20 'C

    sound_speed[1] = config.Read("gas_sound_speed", 340.0f);
    reference_density[1] = config.Read("gas_density", 1.29f);
    gamma[1] = config.Read("gas_gamma", 1.4f);
    viscosity[1] = config.Read("gas_viscosity", 17.9e-6f);
    return false;
}

void FluidProperties::initiate()
{
    float variation_percent = 0.4f;
    for (int i = 0; i < 2; i++)
    {
        coefficient_p2rho[i] = gamma[i] * pow(sound_speed[i], -2.0f) / reference_density[i];
        coefficient_rho2p[i] = 1.0f / coefficient_p2rho[i];
        gamma_inv[i] = 1.0f / gamma[i];
        reference_density_inverse[i] = 1.0f / reference_density[i];

        // min_pressure[i] = (pow((1.0f - variation_percent), gamma[i]) - 1.0f) * coefficient_rho2p[i];
        min_density[i] = (1.0f - variation_percent) * reference_density[i];
    }

}

void SolidProperties::initiate()
{
    switch (constitution_type)
    {
    case STRESS_2D:
        lambda = youngs_modulus * poisson_ratio / (1.0f - poisson_ratio) / (1.0f + poisson_ratio);
        nu = youngs_modulus / (1.0f + poisson_ratio) / 2.0f;
        auxetic_lambda = auxetic_modulus * auxetic_poisson_ratio / (1.0f - auxetic_poisson_ratio) / (1.0f + auxetic_poisson_ratio);
        auxetic_nu = auxetic_modulus / (1.0f + auxetic_poisson_ratio) / 2.0f;
        break;
    case STRAIN_2D:
        lambda = youngs_modulus * poisson_ratio / (1.0f - 2.0f * poisson_ratio) / (1.0f + poisson_ratio);
        nu = youngs_modulus / (1.0f + poisson_ratio) / 2.0f;
        auxetic_lambda = auxetic_modulus * auxetic_poisson_ratio / (1.0f - 2.0f * auxetic_poisson_ratio) / (1.0f + auxetic_poisson_ratio);
        auxetic_nu = auxetic_modulus / (1.0f + auxetic_poisson_ratio) / 2.0f;
        break;
    }

    bulk_modulus = youngs_modulus / (3.0f * (1.0f - 2.0f * poisson_ratio));
    sound_speed = sqrtf(bulk_modulus / reference_density);
}

bool SolidProperties::readNumber(const std::string &path)
{
    std::ifstream infile;
    infile.open(path + "\\solid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read solid particles parameter" << std::endl;
        return true;
    }
    std::string str;
    std::getline(infile, str);
    std::istringstream string_to_num(str);
    string_to_num >> number;
    if (number == 0)
        return false;
    int lines_check(0);
    while (std::getline(infile, str))
        lines_check++;
    if (lines_check != number)
    {
        std::cout << "mismatch solid particle numbers" << std::endl;
        return true;
    }
    infile.close();
    return false;
}

bool SolidProperties::readParameter(const Config &config)
{
    reference_density = config.Read("solid_density", 2700.0f);
    youngs_modulus = config.Read("youngs_modulus", 67.5e9f);
    poisson_ratio = config.Read("poisson_ratio", 0.34f);
    artificial_viscocity[0] = config.Read("arti_vis_alpha", 2.5f);
    artificial_viscocity[1] = config.Read("arti_vis_beta", 2.5f);
    fixed_speed[0] = config.Read("solid_fixed_speed_x", 0.0);
    fixed_speed[1] = config.Read("solid_fixed_speed_y", 0.0);

    bool constitution_type_bool = false;
    constitution_type_bool = config.Read("is_2D_strain", false);
    if (constitution_type_bool)
        constitution_type = STRAIN_2D;
    else
        constitution_type = STRESS_2D;

    bool plastic_type_bool = false;
    plastic_type_bool = config.Read("is_plastic", false);
    if (plastic_type_bool)
        plastic_type = PLASTIC;
    else
        plastic_type = HYPER_ELASTIC;

    auxetic_density = config.Read("auxetic_density", 636.0f);
    auxetic_modulus = config.Read("auxetic_modulus", 0.0f);
    auxetic_poisson_ratio = config.Read("auxetic_poisson_ratio", -0.2f);

    yield_strain = config.Read("yield_strain", 0.0f);
    yield_stress = config.Read("yield_stress", 2.0e7f);
    harden_strain = config.Read("harden_strain", 0.5f);
    harden_coef = config.Read("harden_coef", 100.0f);
    if (auxetic_modulus == 0.0f && yield_strain == 0.0f)
    {
        std::cout << "failed to read solid information" << std::endl;
        return true;
    }
    if (yield_strain != 0.0f)
        auxetic_modulus = yield_stress / yield_strain;
    else
        yield_strain = yield_stress / auxetic_modulus;

    return false;
}

void VirtualProperties::initiate() {}
bool VirtualProperties::readNumber(const std::string &path)
{
    std::ifstream infile;
    infile.open(path + "\\virt_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read virtual particles parameter" << std::endl;
        return true;
    }
    std::string str;
    std::getline(infile, str);
    std::istringstream string_to_num(str);
    string_to_num >> number;
    if (number == 0)
        return false;
    int lines_check(0);
    while (std::getline(infile, str))
        lines_check++;
    if (lines_check != number)
    {
        std::cout << "mismatch virtual particle numbers" << std::endl;
        return true;
    }
    infile.close();
    return false;
}

bool VirtualProperties::readParameter(const Config &config)
{
    return false;
}
