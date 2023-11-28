#ifndef _HOST_DATA_CUH
#define _HOST_DATA_CUH

#include "parameter.h"
#include "types.h"

#include "Eigen/Dense"
#include "thrust/host_vector.h"

class Fluid_Host_List
{
public:
    thrust::host_vector<Eigen::Vector2f> position, velocity;
    thrust::host_vector<float> density, pressure;
    thrust::host_vector<PhaseType> phase_type;
    thrust::host_vector<int> cell_id;

    thrust::host_vector<FluidBoundaryType> fluid_boundary_type;
    thrust::host_vector<NearParticle> has_near_particle;
    thrust::host_vector<Eigen::Matrix2f> du_dx;
    thrust::host_vector<Eigen::Vector2f> gradient, du_dt;
    thrust::host_vector<float> color, drho_dt;

    void addList(const std::string &str, const Parameter &parameter)
    {
        std::istringstream string_to_num(str);
        int phase_id;
        Eigen::Vector2f temp_pos, temp_velo;
        float temp_dens, temp_pres;
        string_to_num >> temp_pos[0] >> temp_pos[1] >> temp_velo[0] >> temp_velo[1] >> temp_dens >> temp_pres >> phase_id;

        position.push_back(temp_pos);
        velocity.push_back(temp_velo);
        density.push_back(temp_dens);
        pressure.push_back(temp_pres);
        phase_type.push_back(PhaseType(phase_id));
        cell_id.push_back(-1);

        fluid_boundary_type.push_back(FREE_SURFACE);
        has_near_particle.push_back(NO_PARTICLE);
        du_dx.push_back(Eigen::Matrix2f::Zero());
        gradient.push_back(Eigen::Vector2f::Zero());
        du_dt.push_back(Eigen::Vector2f::Zero());
        color.push_back(0.0f);
        drho_dt.push_back(0.0f);
    }
};

class Solid_Host_List
{
public:
    thrust::host_vector<Eigen::Vector2d> position, reference_position, velocity;
    thrust::host_vector<SolidBoundaryType> boundary_type;
    thrust::host_vector<float> pressure;
    thrust::host_vector<Eigen::Matrix2d> max_strain;
    thrust::host_vector<MaterialType> material_type;
    thrust::host_vector<int> cell_id;

    thrust::host_vector<Eigen::Matrix2d> deformation, stress_1st_piola, correction;
    thrust::host_vector<Eigen::Vector2d> du_dt;
    thrust::host_vector<Eigen::Vector2f> normal;
    thrust::host_vector<Eigen::Vector2f> couple_dudt;

    void addList(const std::string &str, const Parameter &parameter)
    {
        std::istringstream string_to_num(str);
        int material_type_id(0), boundary_type_id(0);
        Eigen::Vector2d temp_pos, temp_ref_pos, temp_velo;
        float temp_dens, temp_pres;
        string_to_num >> temp_pos[0] >> temp_pos[1] >> temp_velo[0] >> temp_velo[1] >> temp_ref_pos[0] >> temp_ref_pos[1] >> material_type_id >> boundary_type_id;

        reference_position.push_back(temp_ref_pos);
        position.push_back(temp_pos);
        velocity.push_back(temp_velo);
        boundary_type.push_back(SolidBoundaryType(boundary_type_id));
        pressure.push_back(0.0f);
        max_strain.push_back(-parameter.solid.yield_strain * Eigen::Matrix2d::Identity());
        material_type.push_back(MaterialType(material_type_id));
        cell_id.push_back(-1);

        deformation.push_back(Eigen::Matrix2d::Identity());
        stress_1st_piola.push_back(Eigen::Matrix2d::Zero());
        correction.push_back(Eigen::Matrix2d::Identity());
        du_dt.push_back(Eigen::Vector2d::Zero());
        couple_dudt.push_back(Eigen::Vector2f::Zero());
        normal.push_back(Eigen::Vector2f::Zero());
    }
};

class Virtual_Host_List
{
public:
    thrust::host_vector<Eigen::Vector2f> position, normal, velocity;
    thrust::host_vector<float> pressure, density;
    thrust::host_vector<PhaseType> phase_type;
    thrust::host_vector<BoundaryType> boundary_type;
    thrust::host_vector<int> cell_id;

    void addList(const std::string &str, const Parameter &parameter)
    {
        std::istringstream string_to_num(str);

        int boundary_type_id(0);
        Eigen::Vector2f temp_pos;
        string_to_num >> temp_pos[0] >> temp_pos[1] >> boundary_type_id;

        position.push_back(temp_pos);
        normal.push_back(Eigen::Vector2f::Zero());
        velocity.push_back(Eigen::Vector2f::Zero());
        pressure.push_back(0.0f);
        density.push_back(parameter.fluid.reference_density[0]);
        phase_type.push_back(LIQUID);
        boundary_type.push_back(BoundaryType(boundary_type_id));
        cell_id.push_back(-1);
    }
};

#endif