#include "Eigen/Dense"

#include <vector>
#include <fstream>
#include <iostream>

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

enum MaterialType
{
    LINEAR_ELASTIC,
    AUXETIC
};

enum SolidBoundaryType
{
    UNCONSTRAINED_PARTICLE,
    FIXED_PARTICLE
};

class Particle
{
public:
    Eigen::Vector2d position;
};

class Fluid_Particle : public Particle
{
public:
    Eigen::Vector2d velocity;
    double density, pressure;
    PhaseType phase;
};

class Solid_Particle : public Particle
{
public:
    Eigen::Vector2d velocity;
    Eigen::Vector2d reference_position;
    MaterialType material_type;
    SolidBoundaryType boundary_type;
};

class Virtual_Particle : public Particle
{
public:
    BoundaryType boundary_type;
};

const double Pi = 3.1415926;
const double thick = 0.08;
const double half_L = 0.6;
const double middle = 1.5;
const double alpha = 10.0;

bool isInsideWedge(const Eigen::Vector2d &pos, double &shift)
{
    if (pos(1) < shift)
        return false;

    if (pos(0) < middle - half_L || pos(0) > middle + half_L)
        return false;

    double tan_alpha = tan(alpha * Pi / 180.0);
    double top = tan_alpha * half_L + shift + thick;
    if (pos(1) > top)
        return false;
    if (pos(0) < middle)
    {
        double pos_y = -tan_alpha * (pos(0) - middle) + shift;
        if (pos(1) > pos_y && pos(1) < pos_y + thick)
            return true;
        else
            return false;
    }

    if (pos(0) >= middle)
    {
        double pos_y = tan_alpha * (pos(0) - middle) + shift;
        if (pos(1) > pos_y && pos(1) < pos_y + thick)
            return true;
        else
            return false;
    }
    return false;
}

int main()
{
    double dx = 2.5e-3;
    double xmin(0.0), xmax(3.0), ymin(0.0), ymax(3.5);

    double wedge_height = 0.5;
    double water_level = 2.0;
    double shift = water_level + wedge_height;

    std::vector<Fluid_Particle> fluid_list;
    std::vector<Virtual_Particle> virt_list;
    std::vector<Solid_Particle> solid_list;

    int bound_layer = 8;
    int nx = round((xmax - xmin) / dx);
    int ny = round((ymax - ymin) / dx);

    double ref_fluid_density[2];
    ref_fluid_density[0] = 1000.0;
    ref_fluid_density[1] = 1.27;
    double grav = 9.81;
    double gamma[2];
    gamma[0] = 7.0;
    gamma[1] = 1.4;

    int fluid_count[2];
    fluid_count[0] = 0;
    fluid_count[1] = 0;

    for (int i = -bound_layer; i < nx + bound_layer; i++)
    {
        for (int j = -bound_layer; j < ny + bound_layer; j++)
        {
            Eigen::Vector2d pos = Eigen::Vector2d(((double)i + 0.5) * dx + xmin, ((double)j + 0.5) * dx + ymin);
            if (pos(0) < xmin || pos(0) > xmax || pos(1) < ymin || pos(1) > ymax)
            {
                Virtual_Particle virt;
                virt.position = pos;
                virt.boundary_type = SLIP;
                virt_list.push_back(virt);
                continue;
            }

            if (isInsideWedge(pos, shift))
                continue;

            if (pos(1) < water_level)
            {
                Fluid_Particle fluid;
                fluid.position = pos;
                fluid.velocity = Eigen::Vector2d::Zero();
                fluid.pressure = 0.0;
                fluid.density = ref_fluid_density[0];
                fluid.phase = LIQUID;
                fluid_count[0]++;
                fluid_list.push_back(fluid);
            }
            else
            {
                Fluid_Particle fluid;
                fluid.position = pos;
                fluid.velocity = Eigen::Vector2d::Zero();
                fluid.pressure = 0.0;
                fluid.density = ref_fluid_density[1];
                fluid.phase = GAS;
                fluid_count[1]++;
                fluid_list.push_back(fluid);
            }
        }
    }

    double tan_alpha = tan(alpha * Pi / 180.0);
    double velo = -30.0;
    int nx_solid = round((2.0 * half_L) / dx);
    int ny_solid = round((thick) / dx);
    for (int i = 0; i < nx_solid; i++)
    {
        for (int j = 0; j < ny_solid; j++)
        {
            double pos_x = ((double)i + 0.5) * dx + middle - half_L;
            double pos_y;
            if (pos_x < middle)
                pos_y = -tan_alpha * (pos_x - middle) + shift + ((double)j + 0.5) * dx;
            else
                pos_y = tan_alpha * (pos_x - middle) + shift + ((double)j + 0.5) * dx;
            Eigen::Vector2d pos = Eigen::Vector2d(pos_x, pos_y);
            Solid_Particle solid;
            solid.position = pos;
            solid.reference_position = pos;
            if (j * dx < 0.04)
                solid.material_type = AUXETIC;
            else
                solid.material_type = LINEAR_ELASTIC;
            solid.velocity = Eigen::Vector2d(0.0, velo);
            solid.boundary_type = UNCONSTRAINED_PARTICLE;
            if (solid.material_type == LINEAR_ELASTIC)
            {
                if (pos_x < middle - half_L + dx || pos_x > middle + half_L - dx)
                {
                    solid.boundary_type = FIXED_PARTICLE;
                }
                if (pos_x < middle + dx && pos_x > middle - dx && j == ny_solid - 1)
                {
                    solid.boundary_type = FIXED_PARTICLE;
                }
            }
            solid_list.push_back(solid);
        }
    }
    std::fstream fluid_file;
    fluid_file.precision(15);
    fluid_file.open("..\\data\\fluid_particle.dat", std::ios::out);
    fluid_file << fluid_count[0] << " " << fluid_count[1] << std::endl;
    for (auto fp : fluid_list)
    {
        fluid_file << fp.position[0] << " ";
        fluid_file << fp.position[1] << " ";
        fluid_file << fp.velocity[0] << " ";
        fluid_file << fp.velocity[1] << " ";
        fluid_file << fp.density << " ";
        fluid_file << fp.pressure << " ";
        fluid_file << fp.phase << std::endl;
    }
    fluid_file.close();

    std::fstream solid_file;
    solid_file.open("..\\data\\solid_particle.dat", std::ios::out);
    solid_file << solid_list.size() << std::endl;
    for (auto sp : solid_list)
    {
        solid_file << sp.position[0] << " ";
        solid_file << sp.position[1] << " ";
        solid_file << sp.velocity[0] << " ";
        solid_file << sp.velocity[1] << " ";
        solid_file << sp.reference_position[0] << " ";
        solid_file << sp.reference_position[1] << " ";
        solid_file << sp.material_type << " ";
        solid_file << sp.boundary_type << std::endl;
    }
    solid_file.close();

    std::fstream virt_file;
    virt_file.open("..\\data\\virt_particle.dat", std::ios::out);
    virt_file << virt_list.size() << std::endl;
    for (auto vp : virt_list)
    {
        virt_file << vp.position[0] << " ";
        virt_file << vp.position[1] << " ";
        virt_file << vp.boundary_type << std::endl;
    }
    virt_file.close();

    return 0;
}