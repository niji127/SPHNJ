#ifndef _PROBE_H
#define _PROBE_H

#include "sim.cuh"

#include "Eigen/Dense"

class Probe
{
public:
    int id;
    Eigen::Vector2f position, velocity;
    float pressure, density;
    ParticleType particle_type;
    PhaseType phase_type;
    float other_output[4];
    Probe()
    {
        id = 0;
        velocity.setZero();
        pressure = 0.0f;
        particle_type = FLUID;
        phase_type = LIQUID;
        for (int i = 0; i < 4; i++)
            other_output[i] = 0.0f;
    }
    Probe(const Eigen::Vector2f &pos, const int &probe_id) : position(pos), id(probe_id)
    {
        velocity.setZero();
        pressure = 0.0f;
        particle_type = FLUID;
        phase_type = LIQUID;
        for (int i = 0; i < 4; i++)
            other_output[i] = 0.0f;
    }
    void getData();
    void outputData();
};

#endif