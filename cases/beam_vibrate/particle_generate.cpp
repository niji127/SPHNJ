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

double basisFunction(const double &pos_x, const double &length)
{
    double kwL = 1.875;
    double kw = kwL / length;
    return (cos(kwL) + cosh(kwL)) * (cosh(kw * pos_x) - cos(kw * pos_x)) + (sin(kwL) - sinh(kwL)) * (sinh(kw * pos_x) - sin(kw * pos_x));
}

double calculateVelocity(const double &pos_x, const double &length)
{
    double ksi = 0.01f;
    double bulk = 3.25e6;

    bulk = 2.0e6 / (3.0 * (1.0 - 2.0 * 0.4));
    double shear = 7.15e5;
    double rho = 1000.0f;
    double cs = sqrtf(bulk / rho);
    double fL = basisFunction(length, length);
    double fx = basisFunction(pos_x, length);
    return ksi * cs * fx / fL;
}

int main()
{
    double dx = 1.0e-3;
    double xmin(0.0), xmax(0.24), ymin(-0.012), ymax(0.012);

    std::vector<Solid_Particle> solid_list;

    int bound_layer = 4;
    int nx = round((xmax - xmin) / dx);
    int ny = round((ymax - ymin) / dx);

    for (int i = -bound_layer; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            Eigen::Vector2d pos = Eigen::Vector2d(((double)i + 0.5) * dx + xmin, ((double)j + 0.5) * dx + ymin);
            Solid_Particle solid;
            solid.reference_position = pos;
            solid.position = pos;
            solid.material_type = LINEAR_ELASTIC;
            solid.velocity = {0.0f, calculateVelocity(pos(0), xmax - xmin)};

            solid_list.push_back(solid);
        }
    }

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

    std::fstream fluid_file;
    fluid_file.open("..\\data\\fluid_particle.dat", std::ios::out);
    fluid_file << "0 0" << std::endl;
    fluid_file.close();

    std::fstream virt_file;
    virt_file.open("..\\data\\virt_particle.dat", std::ios::out);
    virt_file << 0 << std::endl;
    virt_file.close();

    return 0;
}