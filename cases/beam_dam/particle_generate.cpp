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
    AUXETIC_UNCRUSHED,
    AUXETIC_CRUSHED
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
};

class Virtual_Particle : public Particle
{
public:
    BoundaryType boundary_type;
};

int main()
{
    double dx = 0.00125;
    double xmin(-0.3), xmax(0.1), ymin(0.0), ymax(0.2);

    double solid_left(-0.005);
    double solid_right(0.0);
    double solid_up(0.15);
    double solid_down(0.0);
    double water_level(0.14);

    std::vector<Fluid_Particle> fluid_list;
    std::vector<Virtual_Particle> virt_list;
    std::vector<Solid_Particle> solid_list;

    int bound_layer = 8;
    int nx = round((xmax - xmin) / dx);
    int ny = round((ymax - ymin) / dx);

    double ref_fluid_density[2];
    ref_fluid_density[0] = 1000.0;
    ref_fluid_density[1] = 1.0;
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
            if (pos(0) < solid_left || pos(0) > solid_right || pos(1) < solid_down || pos(1) > solid_up)
            {
                if (pos(0) < xmin || pos(0) > xmax || pos(1) < ymin || pos(1) > ymax)
                {
                    Virtual_Particle virt;
                    virt.position = pos;
                    virt.boundary_type = SLIP;
                    virt_list.push_back(virt);
                }
                else
                {
                    Fluid_Particle fluid;
                    fluid.position = pos;
                    fluid.velocity = Eigen::Vector2d::Zero();
                    fluid.pressure = 0.0;
                    if (pos(0) > solid_right && pos(1) < water_level)
                    {
                        fluid.density = ref_fluid_density[0];
                        fluid.phase = LIQUID;
                        fluid_count[0]++;
                    }
                    else
                    {
                        fluid.density = ref_fluid_density[1];
                        fluid.phase = GAS;
                        fluid_count[1]++;
                    }
                    fluid_list.push_back(fluid);
                }
            }
        }
    }

    int nx_solid = round((solid_right - solid_left) / dx);
    int ny_solid = round((solid_up - solid_down) / dx);
    for (int i = 0; i < nx_solid;i++)
    {
        for (int j = 0; j < ny_solid; j++)
        {
            Eigen::Vector2d pos = Eigen::Vector2d(((double)i + 0.5) * dx + solid_left, ((double)j + 0.5) * dx + solid_down);
            Solid_Particle solid;
            solid.position = pos;
            solid.reference_position = pos;
            solid.material_type = LINEAR_ELASTIC;
            solid.velocity = Eigen::Vector2d::Zero();
            solid_list.push_back(solid);
        }
    }
    std::fstream fluid_file;
    fluid_file.precision(15);
    fluid_file.open("..\\data\\fluid_particle.dat", std::ios::out);
    fluid_file << fluid_count[0] <<" "<<fluid_count[1] <<std::endl;
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
        solid_file << sp.material_type << std::endl;
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