#include "sim.cuh"
#include "convert.h"

#include "host_data.cuh"
#include "thrust/host_vector.h"

#include <string>

Parameter Sim::parameter;
Data_Pointer Sim::data_pointer;

void Sim::showMemoryUsed()
{
    // particles
    int fluid_size = fluid_list.getMemoryUsed();
    int solid_size = solid_list.getMemoryUsed();
    int virt_size = virt_list.getMemoryUsed();

    // cells
    thrust::host_vector<Cell> host_cell;
    host_cell.push_back(cell_list[0]);
    int cell_size = host_cell[0].getMemoryUsed();
    cell_size *= cell_list.size();

    int parameter_size = sizeof(Parameter);
    int sim_size = sizeof(Sim);
    int miscellaneous_size = parameter_size + sim_size;

    int all_size = fluid_size + solid_size + virt_size + cell_size + miscellaneous_size;

    int digit = 2;
    std::cout << "\nAll memory: " << convertUnit(all_size, digit) << std::endl;
    std::cout << "fluid particle: " << convertUnit(fluid_size, digit) << std::endl;
    std::cout << "solid particle: " << convertUnit(solid_size, digit) << std::endl;
    std::cout << "virtual particle: " << convertUnit(virt_size, digit) << std::endl;
    std::cout << "block: " << convertUnit(cell_size, digit) << std::endl;
    std::cout << "other: " << convertUnit(miscellaneous_size, digit) << std::endl;
}

void Sim::readParticleData()
{
    std::cout << "\nReading Particles..." << std::endl;
    bool error_flag(false);
    std::string path = "..\\..\\data";

    clock_t start, end;
    start = clock();
    fluid_list.initiate(parameter);
    solid_list.initiate(parameter);
    virt_list.initiate(parameter);

    error_flag = readFluidData(path);
    error_flag = readSolidData(path) || error_flag;
    error_flag = readVirtData(path) || error_flag;
    end = clock();
    printf("time=%dms\n", end - start);

    if (error_flag)
        throw std::invalid_argument("FAILED TO READ PARTICLES DATA");
}

bool Sim::readFluidData(const std::string &path)
{
    int fluid_num = parameter.fluid.number_total;
    if (fluid_num == 0)
    {
        std::cout << "\nno fluid particle" << std::endl;
        return false;
    }

    if (fluid_num < 0)
    {
        std::cout << "fluid particle number should >= 0" << std::endl;
        return true;
    }

    std::ifstream infile;
    infile.open(path + "\\fluid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read fluid particles data" << std::endl;
        return true;
    }

    Fluid_Host_List fluid_host_list;
    std::string str;
    std::getline(infile, str);
    while (std::getline(infile, str))
        fluid_host_list.addList(str, parameter);
    infile.close();

    fluid_list.copyData(fluid_host_list);

    return false;
}

bool Sim::readSolidData(const std::string &path)
{
    int solid_num = parameter.solid.number;
    if (solid_num == 0)
    {
        std::cout << "\nno solid particle" << std::endl;
        return false;
    }
    if (solid_num < 0)
    {
        std::cout << "solid particle number should >= 0" << std::endl;
        return true;
    }
    std::ifstream infile;
    infile.open(path + "\\solid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read solid particles data" << std::endl;
        return true;
    }

    Solid_Host_List solid_host_list;
    std::string str;
    std::getline(infile, str);
    while (std::getline(infile, str))
        solid_host_list.addList(str, parameter);
    infile.close();

    solid_list.copyData(solid_host_list);

    return false;
}

bool Sim::readVirtData(const std::string &path)
{
    int virt_num = parameter.virt.number;
    if (virt_num == 0)
    {
        std::cout << "\nno virtual particle" << std::endl;
        return false;
    }

    if (virt_num < 0)
    {
        std::cout << "virtual particle number should >= 0" << std::endl;
        return true;
    }

    std::ifstream infile;
    infile.open(path + "\\virt_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read virtual particles data" << std::endl;
        return true;
    }
    Virtual_Host_List virt_host_list;
    std::string str;
    std::getline(infile, str);
    while (std::getline(infile, str))
        virt_host_list.addList(str, parameter);
    infile.close();
    virt_list.copyData(virt_host_list);

    return false;
}

void Sim::createCell()
{
    int cell_number = parameter.domain.cell_number_total;
    int(&cell_n)[2] = parameter.domain.cell_number;
    float(&dx)[2] = parameter.domain.interval;

    
    for (int i = 0; i < cell_n[0]; i++)
    {
        for (int j = 0; j < cell_n[1]; j++)
        {
            Cell temp_cell;
            int id = i * cell_n[1] + j;
            temp_cell.neighbor_cell_id[4] = id;
            if (i != 0)
                temp_cell.neighbor_cell_id[1] = id - cell_n[1];
            if (i != cell_n[0] - 1)
                temp_cell.neighbor_cell_id[7] = id + cell_n[1];
            if (j != 0)
                temp_cell.neighbor_cell_id[3] = id - 1;
            if (j != cell_n[1] - 1)
                temp_cell.neighbor_cell_id[5] = id + 1;

            if ((i != 0) && (j != 0))
                temp_cell.neighbor_cell_id[0] = id - 1 - cell_n[1];
            if ((i != 0) && (j != cell_n[1] - 1))
                temp_cell.neighbor_cell_id[2] = id + 1 - cell_n[1];
            if ((i != cell_n[0] - 1) && (j != 0))
                temp_cell.neighbor_cell_id[6] = id - 1 + cell_n[1];
            if ((i != cell_n[0] - 1) && (j != cell_n[1] - 1))
                temp_cell.neighbor_cell_id[8] = id + 1 + cell_n[1];

            cell_list.push_back(temp_cell);
        }
    }
}