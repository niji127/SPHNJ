#ifndef _PROPERTIES_H
#define _PROPERTIES_H

#include "config.h"
#include "types.h"

class FluidProperties
{
public:
    int number[2];
    int number_total;
    float sound_speed[2], reference_density[2];
    float reference_density_inverse[2];
    float gamma[2], viscosity[2]; // dynamic viscosity
    float particle_shift_coef;
    float gamma_inv[2];
    float coefficient_p2rho[2];
    float coefficient_rho2p[2];
    float min_pressure;
    float min_density[2];
    bool is_fluid_viscid;

    bool readNumber(const std::string &path);
    bool readParameter(const Config &config);
    void initiate();
};

class SolidProperties
{
public:
    int number;
    float reference_density, youngs_modulus, poisson_ratio;
    float lambda, nu, bulk_modulus;
    float artificial_viscocity[2],sound_speed;
    double fixed_speed[2];
    ConstitutionType constitution_type;
    PlasticType plastic_type;

    float auxetic_density;
    float auxetic_modulus;
    float auxetic_poisson_ratio;

    float yield_strain;
    float yield_stress;
    float harden_strain;
    float harden_coef;

    float auxetic_lambda;
    float auxetic_nu;

    bool readNumber(const std::string &path);
    bool readParameter(const Config &config);
    void initiate();
};

class VirtualProperties
{
public:
    int number;

    bool readNumber(const std::string &path);
    bool readParameter(const Config &config);
    void initiate();
};

#endif