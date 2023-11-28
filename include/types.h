#ifndef _TYPES_H
#define _TYPES_H

enum ParticleType
{
    FLUID,
    SOLID,
    VIRT
};

enum PhaseType
{
    LIQUID,
    GAS
};

enum BoundaryType
{
    SLIP,
    NO_SLIP,
    MOVING
};

enum FluidBoundaryType
{
    SPLASH,
    FREE_SURFACE,
    VICINITY,
    INNER
};

enum NearParticle
{
    NO_PARTICLE,
    HAS_PARTICLE,
    HAS_NEAR_PARTICLE,
};

enum MaterialType
{
    LINEAR_ELASTIC,
    AUXETIC
};

enum PlasticType
{
    HYPER_ELASTIC,
    PLASTIC
};

enum SolidBoundaryType
{
    UNCONSTRAINED_PARTICLE,
    FIXED_PARTICLE
};

enum ConstitutionType
{
    STRESS_2D,
    STRAIN_2D
};

enum TransientType
{
    TRANSIENT,
    STEADY
};

struct ParticleID
{
    int particle_id;
    ParticleType particle_type;
};

#endif