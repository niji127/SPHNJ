#include "output.cuh"

#include "thrust/host_vector.h"

#include <windows.h>

void Output::writeFluidVelocity(Sim *sim, const std::string &path, SimTime *time)
{
    int fluid_num = Sim::parameter.fluid.number_total;
    thrust::host_vector<Eigen::Vector2f> fluid_position = sim->fluid_list.position;
    thrust::host_vector<Eigen::Vector2f> fluid_velocity = sim->fluid_list.velocity;

    std::fstream output;
    output.open(path + "\\fluid_velocity-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "Velocity vector file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << fluid_num << " float" << std::endl;
    for (int i = 0; i < fluid_num; i++)
        output << fluid_position[i](0) << " " << fluid_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << fluid_num << std::endl;
    output << "VECTORS velocity float " << std::endl;
    for (int i = 0; i < fluid_num; i++)
        output << fluid_velocity[i](0) << " " << fluid_velocity[i](1) << " 0.0\n";
    output.close();
}

void Output::writeFluidDensity(Sim *sim, const std::string &path, SimTime *time)
{
    int fluid_num = Sim::parameter.fluid.number_total;
    thrust::host_vector<Eigen::Vector2f> fluid_position = sim->fluid_list.position;
    thrust::host_vector<float> fluid_density = sim->fluid_list.density;

    thrust::host_vector<FluidBoundaryType> fluid_type = sim->fluid_list.fluid_boundary_type;
    thrust::host_vector<PhaseType> fluid_phase = sim->fluid_list.phase_type;

    std::fstream output;
    output.open(path + "\\fluid_density-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "scalar file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << fluid_num << " float" << std::endl;
    for (int i = 0; i < fluid_num; i++)
        output << fluid_position[i](0) << " " << fluid_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << fluid_num << std::endl;
    output << "SCALARS density float 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < fluid_num; i++)
    {
        if (fluid_phase[i] == LIQUID)
        {
            // output << fluid_density[i] << "\n";
            output << fluid_type[i] << "\n";
        }
        else
        {
            output << -1 << "\n";
        }
    }
    output.close();
}

void Output::writeFluidPressure(Sim *sim, const std::string &path, SimTime *time)
{
    int fluid_num = Sim::parameter.fluid.number_total;
    thrust::host_vector<Eigen::Vector2f> fluid_position = sim->fluid_list.position;
    thrust::host_vector<float> fluid_pressure = sim->fluid_list.pressure;

    std::fstream output;
    output.open(path + "\\fluid_pressure-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "scalar file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << fluid_num << " float" << std::endl;
    for (int i = 0; i < fluid_num; i++)
        output << fluid_position[i](0) << " " << fluid_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << fluid_num << std::endl;
    output << "SCALARS pressure float 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < fluid_num; i++)
        output << fluid_pressure[i] << "\n";
    output.close();
}

void Output::writeSolidVelocity(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    thrust::host_vector<Eigen::Vector2d> solid_position = sim->solid_list.position;
    thrust::host_vector<Eigen::Vector2d> solid_velocity = sim->solid_list.velocity;
    thrust::host_vector<Eigen::Matrix2d> solid_stress = sim->solid_list.stress_1st_piola;
    thrust::host_vector<Eigen::Matrix2d> solid_max_strain = sim->solid_list.max_strain;

    thrust::host_vector<Eigen::Matrix2d> solid_deform = sim->solid_list.deformation;

    std::fstream output;
    output.open(path + "\\solid_velocity-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "Velocity vector file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << solid_num << " float" << std::endl;
    for (int i = 0; i < solid_num; i++)
        output << solid_position[i](0) << " " << solid_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << solid_num << std::endl;
    output << "VECTORS velocity float " << std::endl;
    for (int i = 0; i < solid_num; i++)
        output << (float)(solid_stress[i](1, 1)) << " " << (float)(solid_max_strain[i](0, 0)) << " " << (float)(solid_max_strain[i](1, 1)) << "\n";
    output.close();
}

void Output::writeSolidPressure(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    thrust::host_vector<Eigen::Vector2d> solid_position = sim->solid_list.position;
    thrust::host_vector<float> solid_pressure = sim->solid_list.pressure;

    std::fstream output;
    output.open(path + "\\solid_pressure-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "scalar file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << solid_num << " float" << std::endl;
    for (int i = 0; i < solid_num; i++)
        output << solid_position[i](0) << " " << solid_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << solid_num << std::endl;
    output << "SCALARS pressure float 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < solid_num; i++)
        output << solid_pressure[i] << "\n";
    output.close();
}

void Output::writeSolidDensity(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    thrust::host_vector<Eigen::Vector2d> solid_position = sim->solid_list.position;
    thrust::host_vector<MaterialType> solid_material = sim->solid_list.material_type;
    float linear_density = Sim::parameter.solid.reference_density;
    float auxetic_density = Sim::parameter.solid.auxetic_density;

    std::fstream output;
    output.open(path + "\\solid_density-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "scalar file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << solid_num << " float" << std::endl;
    for (int i = 0; i < solid_num; i++)
        output << solid_position[i](0) << " " << solid_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << solid_num << std::endl;
    output << "SCALARS density float 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < solid_num; i++)
    {
        if (solid_material[i] == LINEAR_ELASTIC)
            output << linear_density << "\n";
        else
            output << auxetic_density << "\n";
    }
    output.close();
}

void Output::writeSolidStress(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    thrust::host_vector<Eigen::Vector2d> solid_position = sim->solid_list.position;
    thrust::host_vector<Eigen::Matrix2d> solid_stress = sim->solid_list.stress_1st_piola;
    thrust::host_vector<MaterialType> solid_material_type = sim->solid_list.material_type;

    std::fstream output;
    output.open(path + "\\solid_stress-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "Velocity vector file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << solid_num << " float" << std::endl;
    for (int i = 0; i < solid_num; i++)
        output << solid_position[i](0) << " " << solid_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << solid_num << std::endl;
    output << "VECTORS stress float " << std::endl;
    for (int i = 0; i < solid_num; i++)
    {
        if (solid_material_type[i] == LINEAR_ELASTIC)
            output << "0.0 0.0 0.0\n";
        else
            output << (float)(solid_stress[i](0, 0)) << " " << (float)(solid_stress[i](1, 1)) << " " << (float)(solid_stress[i](0, 1)) << "\n";
    }
    output.close();
}

void Output::writeVirtNormal(Sim *sim, const std::string &path, SimTime *time)
{
    int virt_num = Sim::parameter.virt.number;
    thrust::host_vector<Eigen::Vector2f> virt_position = sim->virt_list.position;
    thrust::host_vector<Eigen::Vector2f> virt_normal = sim->virt_list.normal;

    std::fstream output;
    output.open(path + "\\virt_normal-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "Velocity vector file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << virt_num << " float" << std::endl;
    for (int i = 0; i < virt_num; i++)
        output << virt_position[i](0) << " " << virt_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << virt_num << std::endl;
    output << "VECTORS normal float " << std::endl;
    for (int i = 0; i < virt_num; i++)
        output << virt_normal[i](0) << " " << virt_normal[i](1) << " 0.0\n";
    output.close();
}

void Output::writeVirtPressure(Sim *sim, const std::string &path, SimTime *time)
{
    int virt_num = Sim::parameter.virt.number;
    thrust::host_vector<Eigen::Vector2f> virt_position = sim->virt_list.position;
    thrust::host_vector<PhaseType> virt_phase = sim->virt_list.phase_type;

    std::fstream output;
    output.open(path + "\\virt_pressure-" + std::to_string(time->file_step) + ".vtk", std::ios::out);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "scalar file" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << virt_num << " float" << std::endl;
    for (int i = 0; i < virt_num; i++)
        output << virt_position[i](0) << " " << virt_position[i](1) << " 0.0\n";
    output << "POINT_DATA " << virt_num << std::endl;
    output << "SCALARS pressure float 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < virt_num; i++)
        output << virt_phase[i] << "\n";
    output.close();
}

void Output::outputData(Sim *sim)
{
    SimTime *time = Sim::parameter.getTime();
    if (!(time->isOutputData()))
        return;
    std::cout << "Outputing Data..." << std::endl;
    std::cout << "--time=" << time->current_time - time->dt << "s step=" << time->i - 1 << "--" << std::endl;
    std::string file_path = "..\\..\\out";
    if (Sim::parameter.hasFluid())
    {
        std::string fluid_file_path = file_path + "\\fluid";
        CreateDirectory(fluid_file_path.c_str(), NULL);
        writeFluidVelocity(sim, fluid_file_path, time);
        writeFluidDensity(sim, fluid_file_path, time);
        writeFluidPressure(sim, fluid_file_path, time);
    }
    if (Sim::parameter.hasSolid())
    {
        std::string solid_file_path = file_path + "\\solid";
        CreateDirectory(solid_file_path.c_str(), NULL);
        writeSolidVelocity(sim, solid_file_path, time);
        writeSolidPressure(sim, solid_file_path, time);
        writeSolidDensity(sim, solid_file_path, time);
        writeSolidStress(sim, solid_file_path, time);
    }
    if (Sim::parameter.hasVirtual())
    {
        std::string virt_file_path = file_path + "\\virtual";
        CreateDirectory(virt_file_path.c_str(), NULL);
        writeVirtNormal(sim, virt_file_path, time);
        writeVirtPressure(sim, virt_file_path, time);
    }
}

void Output::outputTotalForce(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2f> host_couple_force = sim->solid_list.couple_dudt;
    thrust::host_vector<MaterialType> host_material = sim->solid_list.material_type;
    Eigen::Vector2f gravity = Eigen::Vector2f(Sim::parameter.physics.gravity[0], Sim::parameter.physics.gravity[1]);

    Eigen::Vector2f total_force = Eigen::Vector2f::Zero();

    float rho_elastic = Sim::parameter.solid.reference_density;
    float rho_auxetic = Sim::parameter.solid.auxetic_density;
    float dv = Sim::parameter.kernel.diameter_square;
    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        if (host_material[i] == LINEAR_ELASTIC)
            total_force += host_couple_force[i] * dv * rho_elastic;
        else
            total_force += host_couple_force[i] * dv * rho_auxetic;
    }
    std::fstream output;
    output.open(path + "\\total_force.dat", std::ios::app);
    output << total_force(0) << " " << total_force(1) << std::endl;
    output.close();
}

void Output::outputWedgeDisp(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2d> host_ref_pos = sim->solid_list.reference_position;
    thrust::host_vector<Eigen::Vector2d> host_pos = sim->solid_list.position;

    SimTime *time = Sim::parameter.getTime();
    float current_time = time->current_time - time->dt;

    Eigen::Vector2d current_pos;

    float dx = Sim::parameter.kernel.particle_diameter;

    double middle = 1.5;
    double half_L = 0.6;
    double left_mid_x = middle - 0.5 * half_L;
    double tan_10 = tan(10.0 * PI / 180.0);
    double shift = 2.02;
    double left_mid_y = -tan_10 * (left_mid_x - middle) + shift;
    Eigen::Vector2d middle_wedge(left_mid_x - 0.5 * dx, left_mid_y + 0.5 * dx);
    Eigen::Vector2d down(0.0, -30.0 * current_time);

    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        if ((host_ref_pos[i] - middle_wedge).norm() < 0.5 * dx)
            current_pos = host_pos[i] - host_ref_pos[i] - down;
    }
    std::fstream output;
    output.open(path + "\\wedge_disp.dat", std::ios::app);
    output << current_pos(0) << " " << current_pos(1) << std::endl;
    output.close();
}

void Output::outputEdge(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2f> host_pos = sim->fluid_list.position;
    thrust::host_vector<PhaseType> host_phase = sim->fluid_list.phase_type;

    float frontier = 0.0f;
    for (int i = 0; i < Sim::parameter.fluid.number_total; i++)
    {
        if (host_pos[i](0) > frontier && host_phase[i] == LIQUID)
            frontier = host_pos[i](0);
    }
    std::fstream output;
    output.open(path + "\\frontier.dat", std::ios::app);
    output << frontier + Sim::parameter.kernel.particle_diameter * 0.5f << std::endl;
    output.close();

    float dx = Sim::parameter.kernel.particle_diameter;
    float height_0(-0.5f * dx), height_1(-0.5f * dx);
    float x0(5.366f - 0.825f), x1(5.366f - 1.653f);
    for (int i = 0; i < Sim::parameter.fluid.number_total; i++)
    {
        if ((host_pos[i](1) > height_0) && (host_pos[i](0) < x0 + 0.5f * dx) && (host_pos[i](0) > x0 - 0.5f * dx) && (host_phase[i] == LIQUID))
            height_0 = host_pos[i](1);

        if ((host_pos[i](1) > height_1) && (host_pos[i](0) < x1 + 0.5f * dx) && (host_pos[i](0) > x1 - 0.5f * dx) && (host_phase[i] == LIQUID))
            height_1 = host_pos[i](1);
    }
    output.open(path + "\\height.dat", std::ios::app);
    output << height_0 + dx * 0.5f << " " << height_1 + dx * 0.5f << std::endl;
    output.close();
}

void Output::outputMiddlePosition(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2d> host_ref_pos = sim->solid_list.reference_position;
    thrust::host_vector<Eigen::Vector2d> host_pos = sim->solid_list.position;

    float dx = Sim::parameter.kernel.particle_diameter;
    float xmin = Sim::parameter.domain.domain_xmin[0];
    float xmax = Sim::parameter.domain.domain_xmax[0];
    float ymin = Sim::parameter.domain.domain_xmin[1];
    float ymax = Sim::parameter.domain.domain_xmax[1];
    float thick = 0.05f;
    Eigen::Vector2d end_ref_pos = {0.5 * (xmax + xmin) - 0.5 * dx, -thick + 0.5 * dx};
    Eigen::Vector2d end_pos;
    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        if ((host_ref_pos[i] - end_ref_pos).norm() < 0.5f * dx)
        {
            end_pos = host_pos[i];
            break;
        }
    }
    end_pos += Eigen::Vector2d(0.5 * dx, thick - 0.5 * dx);
    std::fstream output;
    output.open(path + "\\middle_position.dat", std::ios::app);
    output << end_pos(0) << " " << end_pos(1) << std::endl;
    output.close();
}

void Output::outputEndPosition(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2d> host_ref_pos = sim->solid_list.reference_position;
    thrust::host_vector<Eigen::Vector2d> host_pos = sim->solid_list.position;

    float dx = Sim::parameter.kernel.particle_diameter;
    float xmin = Sim::parameter.domain.domain_xmin[0];
    float xmax = Sim::parameter.domain.domain_xmax[0];
    float ymin = Sim::parameter.domain.domain_xmin[1];
    float ymax = Sim::parameter.domain.domain_xmax[1];
    Eigen::Vector2d end_ref_pos = {xmax - 0.5 * dx, 0.5 * (ymin + ymax) - 0.5 * dx};
    Eigen::Vector2d end_pos;
    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        if ((host_ref_pos[i] - end_ref_pos).norm() < 0.5f * dx)
        {
            end_pos = host_pos[i];
            break;
        }
    }
    end_pos += Eigen::Vector2d(0.5 * dx, 0.5 * dx);
    std::fstream output;
    output.open(path + "\\end_position.dat", std::ios::app);
    output << end_pos(0) << " " << end_pos(1) << std::endl;
    output.close();
}

void Output::outputTotalEnergy(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2d> host_ref_pos = sim->solid_list.reference_position;
    thrust::host_vector<Eigen::Vector2d> host_vel = sim->solid_list.velocity;
    thrust::host_vector<Eigen::Matrix2d> host_deform = sim->solid_list.deformation;

    float rho = Sim::parameter.solid.reference_density;
    float lambda = Sim::parameter.solid.lambda;
    float nu = Sim::parameter.solid.nu;
    float dv = Sim::parameter.kernel.diameter_square;

    double Ep(0.0), Ek(0.0);

    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        Eigen::Matrix2d green = 0.5 * (host_deform[i].transpose() * host_deform[i] - Eigen::Matrix2d::Identity());

        Ep += 0.5 * lambda * pow(green.trace(), 2.0) * dv;
        Ep += nu * green.cwiseProduct(green).sum() * dv;

        if (host_ref_pos[i](0) < 0.0)
            continue;

        Ek += 0.5 * rho * dv * host_vel[i].squaredNorm();
    }
    std::fstream output;
    output.open(path + "\\total_energy.dat", std::ios::app);
    output << Ep << " " << Ek << std::endl;
    output.close();
}

void Output::outputTankBeamEnd(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2d> host_ref_pos = sim->solid_list.reference_position;
    thrust::host_vector<Eigen::Vector2d> host_pos = sim->solid_list.position;

    float dx = Sim::parameter.kernel.particle_diameter;
    double thick = 0.005;
    Eigen::Vector2d end_ref_pos = {-0.5 * thick + 0.5 * dx, 0.5 * dx};
    Eigen::Vector2d end_pos;
    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        if ((host_ref_pos[i] - end_ref_pos).norm() < 0.5 * dx)
        {
            end_pos = host_pos[i] - end_ref_pos;
            break;
        }
    }
    std::fstream output;
    output.open(path + "\\beam_tank_end.dat", std::ios::app);
    output << end_pos(0) << " " << end_pos(1) << std::endl;
    output.close();
}

void Output::outputCrushForce(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Vector2d> host_ref_pos = sim->solid_list.reference_position;
    thrust::host_vector<Eigen::Vector2f> host_dudt = sim->solid_list.couple_dudt;

    float dx = Sim::parameter.kernel.particle_diameter;
    float rho = Sim::parameter.solid.auxetic_density;
    float dv = Sim::parameter.kernel.diameter_square;
    float sum = 0.0f;
    float force = 0.0f;
    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        if (host_ref_pos[i](1) < dx && host_ref_pos[i](1) > 0.0)
        {
            force += host_dudt[i](1) * rho * dv;
        }
    }
    std::fstream output;
    output.open(path + "\\crush_force.dat", std::ios::app);
    output << force / 16.2e-3 << std::endl;
    output.close();
}

void Output::outputVonMisesStress(Sim *sim, const std::string &path)
{
    thrust::host_vector<Eigen::Matrix2d> stress_1st_piola = sim->solid_list.stress_1st_piola;
    thrust::host_vector<Eigen::Matrix2d> deformation = sim->solid_list.deformation;
    thrust::host_vector<MaterialType> material_type = sim->solid_list.material_type;

    thrust::host_vector<Eigen::Vector2d> position = sim->solid_list.position;

    double max_stress = 0.0;
    double max_stress_al = 0.0;
    double max_stress_au = 0.0;

    float rho_elastic = Sim::parameter.solid.reference_density;
    float rho_auxetic = Sim::parameter.solid.auxetic_density;
    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        double mises_stress = 0.0;
        Eigen::Matrix2d stress = deformation[i].inverse() * stress_1st_piola[i];
        mises_stress += 0.5 * pow(stress(0, 0) - stress(1, 1), 2.0);
        mises_stress += 3.0 * pow(stress(0, 1), 2.0);
        mises_stress = sqrt(mises_stress);
        if (material_type[i] == LINEAR_ELASTIC)
        {
            if (mises_stress > max_stress_al)
            {
                max_stress_al = mises_stress;
            }
        }
        else if (material_type[i] == AUXETIC)
        {
            if (mises_stress > max_stress_au)
            {
                max_stress_au = mises_stress;
            }
        }
    }
    max_stress = max(max_stress_al, max_stress_au);
    std::fstream output;
    output.open(path + "\\von_mises_stress.dat", std::ios::app);
    output << max_stress_al << " " << max_stress_au << " " << max_stress << std::endl;
    output.close();
}

void Output::outputMeanPressure(Sim *sim, const std::string &path)
{
    thrust::host_vector<float> pressure = sim->solid_list.pressure;

    float total_pressure = 0.0f;
    for (int i = 0; i < Sim::parameter.solid.number; i++)
    {
        total_pressure += pressure[i];
    }
    total_pressure /= (float)Sim::parameter.solid.number;
    std::fstream output;
    output.open(path + "\\mean_pressure.dat", std::ios::app);
    output << total_pressure << std::endl;
    output.close();
}

void Output::outputGlobalData(Sim *sim)
{
    std::string path = "..\\..\\out\\global";
    CreateDirectory(path.c_str(), NULL);
    // outputEdge(sim, path);

    // outputTotalEnergy(sim, path);
    // outputEndPosition(sim, path);

    // outputMiddlePosition(sim, path);

    // outputTankBeamEnd(sim, path);

    outputTotalForce(sim, path);
    outputWedgeDisp(sim, path);

    // outputCrushForce(sim, path);

    outputVonMisesStress(sim, path);
    outputMeanPressure(sim, path);
}