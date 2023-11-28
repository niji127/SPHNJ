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

const double Pi = 3.1415926;
const double thick = 0.018;
const double radius = 0.05;
const double middle = 0.15;
const double dx = 1.0e-3;

bool isInsideCircle(const Eigen::Vector2d &pos, double &shift)
{
    Eigen::Vector2d center(middle, shift);
    if ((pos - center).norm() > radius+dx)
        return false;
    else
        return true;
}

int main()
{
    
    double xmin(0.0), xmax(0.3), ymin(0.0), ymax(0.5);

    double circle_height = 0.1;
    double water_level = 0.3;
    double shift = water_level + circle_height;

    double velo = -50.0;
    double inner_thick = 0.004;
    double middle_thick = 0.01;
    double outer_thick = 0.004;

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

            if (isInsideCircle(pos, shift))
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


    double total_thick = inner_thick + middle_thick + outer_thick;
    int inner_layer = round(inner_thick / dx);
    int middle_layer = round(middle_thick / dx);
    int outer_layer = round(outer_thick / dx);
    Eigen::Vector2d center(middle, shift);
    for (int i = 0; i < inner_layer; i++)
    {
        double r = radius - total_thick + ((double)i + 0.5) * dx;
        int radius_n = round(2.0 * Pi * r / dx);
        double d_theta = 2.0 * Pi * r / (double)radius_n;
        for (int j = 0; j < radius_n; j++)
        {
            double theta = (double)j * d_theta / r;

            double pos_x = r * cos(theta);
            double pos_y = r * sin(theta);
            Eigen::Vector2d pos = Eigen::Vector2d(pos_x, pos_y);
            pos += center;
            Solid_Particle solid;
            solid.position = pos;
            solid.reference_position = pos;
            solid.material_type = LINEAR_ELASTIC;
            solid.velocity = Eigen::Vector2d(0.0, velo);
            solid_list.push_back(solid);
        }
    }
    for (int i = 0; i < middle_layer; i++)
    {
        double r = radius - outer_thick - middle_thick + ((double)i + 0.5) * dx;
        int radius_n = round(2.0 * Pi * r / dx);
        double d_theta = 2.0 * Pi * r / (double)radius_n;
        for (int j = 0; j < radius_n; j++)
        {
            double theta = (double)j * d_theta / r;

            double pos_x = r * cos(theta);
            double pos_y = r * sin(theta);
            Eigen::Vector2d pos = Eigen::Vector2d(pos_x, pos_y);
            pos += center;
            Solid_Particle solid;
            solid.position = pos;
            solid.reference_position = pos;
            solid.material_type = AUXETIC_UNCRUSHED;
            solid.velocity = Eigen::Vector2d(0.0, velo);
            solid_list.push_back(solid);
        }
    }
    for (int i = 0; i < outer_layer; i++)
    {
        double r = radius - outer_thick + ((double)i + 0.5) * dx;
        int radius_n = round(2.0 * Pi * r / dx);
        double d_theta = 2.0 * Pi * r / (double)radius_n;
        for (int j = 0; j < radius_n; j++)
        {
            double theta = (double)j * d_theta / r;

            double pos_x = r * cos(theta);
            double pos_y = r * sin(theta);
            Eigen::Vector2d pos = Eigen::Vector2d(pos_x, pos_y);
            pos += center;
            Solid_Particle solid;
            solid.position = pos;
            solid.reference_position = pos;
            solid.material_type = LINEAR_ELASTIC;
            solid.velocity = Eigen::Vector2d(0.0, velo);
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