#ifndef _KERNEL_H
#define _KERNEL_H

#include "config.h"

class Kernel
{
public:
    float smoothing_length;
    float particle_diameter;
    float particle_volume;
    float diameter_square;

    int thread_num;

    float smoothing_length_square;
    float smoothing_length_inverse;
    float smoothing_length_inverse_square;

    float impact_length_hsml_ratio;

    float impact_distance_square;

    float kernel_coefficient_1;
    float kernel_coefficient_2;

    float kernel_differential;
    float coefficient_3d_to_2d;

    bool readParameter(const Config &config);
    void initiate();
    bool readNumber(const std::string &path);
};

#endif