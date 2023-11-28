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
    Eigen::Vector2f position;
};

class Fluid_Particle : public Particle
{
public:
    Eigen::Vector2f velocity;
    float density, pressure;
    PhaseType phase;
};

class Solid_Particle : public Particle
{
public:
    Eigen::Vector2d velocity, position;
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
    float dx = 5.0e-3;
    float xmin(0.0f), xmax(5.366f), ymin(0.0f), ymax(2.5f);
    Eigen::Vector2f liquid_boundary(2.0f, 1.0f);

    std::vector<Fluid_Particle> fluid_list;
    std::vector<Solid_Particle> solid_list;
    std::vector<Virtual_Particle> virt_list;

    int bound_layer = 8;
    int nx = (int)((xmax - xmin) / dx);
    int ny = (int)((ymax - ymin) / dx);

    double liquid_density = 1.0;
    double gas_density = 1.0e-3;
    int fluid_num[2];
    fluid_num[0] = 0;
    fluid_num[1] = 0;

    for (int i = -bound_layer; i < nx + bound_layer; i++)
    {
        for (int j = -bound_layer; j < ny + bound_layer; j++)
        {
            Eigen::Vector2f pos = Eigen::Vector2f(((float)i + 0.5f) * dx + xmin, ((float)j + 0.5f) * dx + ymin);

            // virtual particles
            if (pos(0) < xmin || pos(0) > xmax || pos(1) < ymin || pos(1) > ymax)
            {
                Virtual_Particle virt;
                virt.position = pos;
                virt.boundary_type = SLIP;
                virt_list.push_back(virt);
            }
            // liquid particles
            else if (pos(0) < liquid_boundary(0) && pos(1) < liquid_boundary(1))
            {
                Fluid_Particle liquid;
                liquid.position = pos;
                liquid.velocity = Eigen::Vector2f::Zero();
                liquid.density = liquid_density;
                liquid.pressure = 0.0f;
                liquid.phase = LIQUID;
                fluid_list.push_back(liquid);
                fluid_num[0]++;
            }
            // gas particles
            else
            {
                Fluid_Particle gas;
                gas.position = pos;
                gas.velocity = Eigen::Vector2f::Zero();
                gas.density = gas_density;
                gas.pressure = 0.0f;
                gas.phase = GAS;
                fluid_list.push_back(gas);
                fluid_num[1]++;
            }
        }
    }

    std::fstream fluid_file;
    fluid_file.open("..\\data\\fluid_particle.dat", std::ios::out);
    fluid_file << fluid_num[0] << " " << fluid_num[1] << std::endl;
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

    std::fstream solid_file;
    solid_file.open("..\\data\\solid_particle.dat", std::ios::out);
    solid_file << 0 << std::endl;
    solid_file.close();

    return 0;
}