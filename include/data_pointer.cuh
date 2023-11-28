#ifndef _DATA_POINTER_CUH
#define _DATA_POINTER_CUH

#include "particle.cuh"
#include "cell.cuh"

class Fluid_Data
{
public:
    Eigen::Vector2f *position, *velocity;
    float *density, *pressure;
    PhaseType *phase_type;
    int *cell_id;

    FluidBoundaryType *fluid_boundary_type;
    NearParticle *has_near_particle;
    Eigen::Matrix2f *du_dx;
    Eigen::Vector2f *gradient, *du_dt;
    float *color, *drho_dt;
};

class Solid_Data
{
public:
    Eigen::Vector2d *position, *reference_position, *velocity;
    SolidBoundaryType *boundary_type;
    float *pressure;
    Eigen::Matrix2d *max_strain;
    MaterialType *material_type;
    int *cell_id;

    Eigen::Matrix2d *deformation,*stress_1st_piola, *correction;
    Eigen::Vector2d *du_dt;
    Eigen::Vector2f *couple_dudt;
    Eigen::Vector2f *normal;
};

class Virtual_Data
{
public:
    Eigen::Vector2f *position, *normal, *velocity;
    float *pressure, *density;
    PhaseType *phase_type;
    BoundaryType *boundary_type;
    int *cell_id;
};

class Data_Pointer
{
public:
    Fluid_Data fluid;
    Solid_Data solid;
    Virtual_Data virt;
    Cell *cell;
    void fluidInitiate(Fluid_List &fluid_list)
    {
        fluid.position = thrust::raw_pointer_cast(fluid_list.position.data());
        fluid.velocity = thrust::raw_pointer_cast(fluid_list.velocity.data());
        fluid.density = thrust::raw_pointer_cast(fluid_list.density.data());
        fluid.pressure = thrust::raw_pointer_cast(fluid_list.pressure.data());
        fluid.phase_type = thrust::raw_pointer_cast(fluid_list.phase_type.data());
        fluid.cell_id = thrust::raw_pointer_cast(fluid_list.cell_id.data());

        fluid.fluid_boundary_type = thrust::raw_pointer_cast(fluid_list.fluid_boundary_type.data());
        fluid.has_near_particle = thrust::raw_pointer_cast(fluid_list.has_near_particle.data());
        fluid.du_dx = thrust::raw_pointer_cast(fluid_list.du_dx.data());
        fluid.gradient = thrust::raw_pointer_cast(fluid_list.gradient.data());
        fluid.du_dt = thrust::raw_pointer_cast(fluid_list.du_dt.data());
        fluid.color = thrust::raw_pointer_cast(fluid_list.color.data());
        fluid.drho_dt = thrust::raw_pointer_cast(fluid_list.drho_dt.data());
    }

    void solidInitiate(Solid_List &solid_list)
    {
        solid.reference_position = thrust::raw_pointer_cast(solid_list.reference_position.data());
        solid.position = thrust::raw_pointer_cast(solid_list.position.data());
        solid.velocity = thrust::raw_pointer_cast(solid_list.velocity.data());
        solid.boundary_type = thrust::raw_pointer_cast(solid_list.boundary_type.data());
        solid.pressure = thrust::raw_pointer_cast(solid_list.pressure.data());
        solid.max_strain = thrust::raw_pointer_cast(solid_list.max_strain.data());
        solid.material_type = thrust::raw_pointer_cast(solid_list.material_type.data());
        solid.cell_id = thrust::raw_pointer_cast(solid_list.cell_id.data());

        solid.deformation = thrust::raw_pointer_cast(solid_list.deformation.data());
        solid.stress_1st_piola = thrust::raw_pointer_cast(solid_list.stress_1st_piola.data());
        solid.correction = thrust::raw_pointer_cast(solid_list.correction.data());
        solid.du_dt = thrust::raw_pointer_cast(solid_list.du_dt.data());
        solid.couple_dudt = thrust::raw_pointer_cast(solid_list.couple_dudt.data());
        solid.normal = thrust::raw_pointer_cast(solid_list.normal.data());
    }

    void virtInitiate(Virtual_List &virt_list)
    {
        virt.position = thrust::raw_pointer_cast(virt_list.position.data());
        virt.normal = thrust::raw_pointer_cast(virt_list.normal.data());
        virt.velocity = thrust::raw_pointer_cast(virt_list.velocity.data());
        virt.pressure = thrust::raw_pointer_cast(virt_list.pressure.data());
        virt.density = thrust::raw_pointer_cast(virt_list.density.data());
        virt.phase_type = thrust::raw_pointer_cast(virt_list.phase_type.data());
        virt.boundary_type = thrust::raw_pointer_cast(virt_list.boundary_type.data());
        virt.cell_id = thrust::raw_pointer_cast(virt_list.cell_id.data());
    }

    void cellInitiate(thrust::device_vector<Cell> &cell_list)
    {
        cell = thrust::raw_pointer_cast(cell_list.data());
    }

    void initiate(Fluid_List &fluid_list, Solid_List &solid_list, Virtual_List &virt_list, thrust::device_vector<Cell> &cell_list)
    {
        fluidInitiate(fluid_list);
        solidInitiate(solid_list);
        virtInitiate(virt_list);
        cellInitiate(cell_list);
    }
};

#endif