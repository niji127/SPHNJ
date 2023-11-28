#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "global.h"
#include "types.h"
#include "parameter.h"
#include "host_data.cuh"

#include "Eigen/Dense"
#include "thrust/device_vector.h"
#include "thrust/sort.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/tuple.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <string>
#include <sstream>

class Fluid_List
{
public:
    thrust::device_vector<Eigen::Vector2f> position, velocity;
    thrust::device_vector<float> density, pressure;
    thrust::device_vector<PhaseType> phase_type;
    thrust::device_vector<int> cell_id;

    thrust::device_vector<FluidBoundaryType> fluid_boundary_type;
    thrust::device_vector<NearParticle> has_near_particle;
    thrust::device_vector<Eigen::Matrix2f> du_dx;
    thrust::device_vector<Eigen::Vector2f> gradient, du_dt;
    thrust::device_vector<float> color, drho_dt;

    thrust::device_vector<int> cell_id_list;
    thrust::device_vector<int> cell_counts;
    thrust::device_vector<int> cell_indices;
    int cell_num;

    size_t getMemoryUsed()
    {
        size_t total_size = 0;
        total_size += position.size() * sizeof(Eigen::Vector2f) + sizeof(position);
        total_size += velocity.size() * sizeof(Eigen::Vector2f) + sizeof(velocity);
        total_size += density.size() * sizeof(float) + sizeof(density);
        total_size += pressure.size() * sizeof(float) + sizeof(pressure);
        total_size += phase_type.size() * sizeof(PhaseType) + sizeof(phase_type);
        total_size += cell_id.size() * sizeof(int) + sizeof(cell_id);

        total_size += fluid_boundary_type.size() * sizeof(FluidBoundaryType) + sizeof(fluid_boundary_type);
        total_size += has_near_particle.size() * sizeof(NearParticle) + sizeof(has_near_particle);
        total_size += du_dx.size() * sizeof(Eigen::Matrix2f) + sizeof(du_dx);
        total_size += gradient.size() * sizeof(Eigen::Vector2f) + sizeof(gradient);
        total_size += du_dt.size() * sizeof(Eigen::Vector2f) + sizeof(du_dt);
        total_size += color.size() * sizeof(float) + sizeof(color);
        total_size += drho_dt.size() * sizeof(float) + sizeof(drho_dt);

        total_size += cell_id_list.size() * sizeof(int) + sizeof(cell_id_list);
        total_size += cell_counts.size() * sizeof(int) + sizeof(cell_counts);
        total_size += cell_indices.size() * sizeof(int) + sizeof(cell_indices);
        total_size += sizeof(cell_num);

        return total_size;
    }

    void copyData(const Fluid_Host_List &host_data)
    {
        position = host_data.position;
        velocity = host_data.velocity;
        density = host_data.density;
        pressure = host_data.pressure;
        phase_type = host_data.phase_type;
        cell_id = host_data.cell_id;

        fluid_boundary_type = host_data.fluid_boundary_type;
        has_near_particle = host_data.has_near_particle;
        du_dx = host_data.du_dx;
        gradient = host_data.gradient;
        du_dt = host_data.du_dt;
        color = host_data.color;
        drho_dt = host_data.drho_dt;
    }

    void initiate(const Parameter &parameter)
    {
        int fluid_num = parameter.fluid.number_total;

        cell_id_list.resize(fluid_num);
        cell_counts.resize(fluid_num);
        cell_indices.resize(fluid_num);
        cell_num = 0;
    }

    auto getZipProperties()
    {
        return thrust::make_zip_iterator(thrust::make_tuple(position.begin(), velocity.begin(), density.begin(), pressure.begin(), phase_type.begin(), fluid_boundary_type.begin()));
    }

    void sort();
};

class Solid_List
{
public:
    thrust::device_vector<Eigen::Vector2d> position, reference_position, velocity;
    thrust::device_vector<SolidBoundaryType> boundary_type;
    thrust::device_vector<Eigen::Matrix2d> max_strain;
    thrust::device_vector<MaterialType> material_type;
    thrust::device_vector<float> pressure;
    thrust::device_vector<int> cell_id;

    thrust::device_vector<Eigen::Matrix2d> deformation, stress_1st_piola, correction;
    thrust::device_vector<Eigen::Vector2d> du_dt;
    thrust::device_vector<Eigen::Vector2f> normal;

    thrust::device_vector<Eigen::Vector2f> couple_dudt;

    thrust::device_vector<int> cell_id_list;
    thrust::device_vector<int> cell_counts;
    thrust::device_vector<int> cell_indices;
    int cell_num;

    size_t getMemoryUsed()
    {
        size_t total_size = 0;
        total_size += reference_position.size() * sizeof(Eigen::Vector2d) + sizeof(reference_position);
        total_size += position.size() * sizeof(Eigen::Vector2d) + sizeof(position);
        total_size += velocity.size() * sizeof(Eigen::Vector2d) + sizeof(velocity);
        total_size += boundary_type.size() * sizeof(SolidBoundaryType) + sizeof(boundary_type);
        total_size += max_strain.size() * sizeof(Eigen::Matrix2d) + sizeof(max_strain);
        total_size += pressure.size() * sizeof(float) + sizeof(pressure);
        total_size += material_type.size() * sizeof(MaterialType) + sizeof(material_type);
        total_size += cell_id.size() * sizeof(int) + sizeof(cell_id);

        total_size += deformation.size() * sizeof(Eigen::Matrix2d) + sizeof(deformation);
        total_size += stress_1st_piola.size() * sizeof(Eigen::Matrix2d) + sizeof(stress_1st_piola);
        total_size += correction.size() * sizeof(Eigen::Matrix2d) + sizeof(correction);
        total_size += du_dt.size() * sizeof(Eigen::Vector2d) + sizeof(du_dt);
        total_size += normal.size() * sizeof(Eigen::Vector2f) + sizeof(normal);
        total_size += couple_dudt.size() * sizeof(Eigen::Vector2f) + sizeof(couple_dudt);

        total_size += cell_id_list.size() * sizeof(int) + sizeof(cell_id_list);
        total_size += cell_counts.size() * sizeof(int) + sizeof(cell_counts);
        total_size += cell_indices.size() * sizeof(int) + sizeof(cell_indices);
        total_size += sizeof(cell_num);

        return total_size;
    }

    void copyData(const Solid_Host_List &host_data)
    {
        position = host_data.position;
        reference_position = host_data.reference_position;
        velocity = host_data.velocity;
        boundary_type = host_data.boundary_type;
        pressure = host_data.pressure;
        max_strain = host_data.max_strain;
        material_type = host_data.material_type;
        cell_id = host_data.cell_id;

        deformation = host_data.deformation;
        stress_1st_piola = host_data.stress_1st_piola;
        correction = host_data.correction;
        du_dt = host_data.du_dt;
        normal = host_data.normal;
        couple_dudt = host_data.couple_dudt;
    }

    void initiate(const Parameter &parameter)
    {
        cell_id_list.resize(parameter.solid.number);
        cell_counts.resize(parameter.solid.number);
        cell_indices.resize(parameter.solid.number);
        cell_num = 0;
    }

    auto getZipProperties()
    {
        return thrust::make_zip_iterator(thrust::make_tuple(reference_position.begin(), position.begin(), velocity.begin(), material_type.begin(), boundary_type.begin()));
    }

    void sort();
};

class Virtual_List
{
public:
    thrust::device_vector<Eigen::Vector2f> position, normal, velocity;
    thrust::device_vector<float> pressure, density;
    thrust::device_vector<PhaseType> phase_type;
    thrust::device_vector<BoundaryType> boundary_type;
    thrust::device_vector<int> cell_id;

    thrust::device_vector<int> cell_id_list;
    thrust::device_vector<int> cell_counts;
    thrust::device_vector<int> cell_indices;
    int cell_num;

    size_t getMemoryUsed()
    {
        size_t total_size = 0;
        total_size += position.size() * sizeof(Eigen::Vector2f) + sizeof(position);
        total_size += normal.size() * sizeof(Eigen::Vector2f) + sizeof(normal);
        total_size += velocity.size() * sizeof(Eigen::Vector2f) + sizeof(velocity);
        total_size += pressure.size() * sizeof(float) + sizeof(pressure);
        total_size += density.size() * sizeof(float) + sizeof(density);
        total_size += phase_type.size() * sizeof(PhaseType) + sizeof(phase_type);
        total_size += boundary_type.size() * sizeof(BoundaryType) + sizeof(boundary_type);
        total_size += cell_id.size() * sizeof(int) + sizeof(cell_id);

        total_size += cell_id_list.size() * sizeof(int) + sizeof(cell_id_list);
        total_size += cell_counts.size() * sizeof(int) + sizeof(cell_counts);
        total_size += cell_indices.size() * sizeof(int) + sizeof(cell_indices);
        total_size += sizeof(cell_num);

        return total_size;
    }

    void copyData(const Virtual_Host_List &host_data)
    {
        position=host_data.position;
        normal = host_data.normal;
        velocity = host_data.velocity;
        pressure = host_data.pressure;
        density = host_data.density;
        phase_type = host_data.phase_type;
        boundary_type = host_data.boundary_type;
        cell_id = host_data.cell_id;
    }

    void initiate(const Parameter &parameter)
    {
        cell_id_list.resize(parameter.virt.number);
        cell_counts.resize(parameter.virt.number);
        cell_indices.resize(parameter.virt.number);
        cell_num = 0;
    }

    auto getZipProperties()
    {
        return thrust::make_zip_iterator(thrust::make_tuple(position.begin(), normal.begin(), boundary_type.begin()));
    }

    void sort();
};
#endif