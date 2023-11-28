#include "parameter.h"

void Parameter::showInfo()
{
    std::cout << "Parameters reading finished" << std::endl;
    std::cout << "\n[physics]" << std::endl;
    std::cout << "gravity=" << physics.gravity[1] << std::endl;
    std::cout << "\n[particles]" << std::endl;
    std::cout << "fluid number:" << fluid.number_total << std::endl;
    std::cout << "liquid:" << fluid.number[0] << "    gas:" << fluid.number[1] << std::endl;
    std::cout << "solid number:" << solid.number << std::endl;
    std::cout << "virtual number:" << virt.number << std::endl;
    std::cout << "\n[domain]" << std::endl;
    std::cout << "x(" << domain.domain_xmin[0] << "," << domain.domain_xmax[0];
    std::cout << ") y(" << domain.domain_xmin[1] << "," << domain.domain_xmax[1] << ")" << std::endl;
    std::cout << "block number:" << domain.cell_number_total << "  ";
    std::cout << domain.cell_number[0] << "*" << domain.cell_number[1] << std::endl;
    std::cout << "\n[time]" << std::endl;
    std::cout << "dt=" << time.dt << "  current step=" << time.i << std::endl;
    std::cout << "start time=" << (float)time.istart * time.dt;
    std::cout << "  end time=" << (float)time.iend * time.dt << std::endl;
}

void Parameter::initiateParameter()
{
    physics.initiate();
    fluid.initiate();
    solid.initiate();
    virt.initiate();
    domain.initiate();
    time.initiate();
    kernel.initiate();
    domain.initiate(kernel);
}

void Parameter::readFieldParameter()
{
    bool error_flag(false);
    std::string path = "..\\..\\config\\input.ini";
    Config config(path);

    error_flag = physics.readParameter(config);
    error_flag = fluid.readParameter(config) || error_flag;
    error_flag = solid.readParameter(config) || error_flag;
    error_flag = virt.readParameter(config) || error_flag;
    error_flag = domain.readParameter(config) || error_flag;
    error_flag = time.readParameter(config) || error_flag;
    error_flag = kernel.readParameter(config) || error_flag;

    if (error_flag)
        throw std::invalid_argument("FAILED TO READ INPUT.INI");
}

void Parameter::readParticleNumber()
{
    bool error_flag(false);
    std::string path = "..\\..\\data";
    error_flag = fluid.readNumber(path);
    error_flag = solid.readNumber(path) || error_flag;
    error_flag = virt.readNumber(path) || error_flag;
    if (error_flag)
        throw std::invalid_argument("FAILED TO READ PARTICLES PARAMETER");
}

bool Parameter::hasFluid()
{
    return fluid.number_total == 0 ? false : true;
}

bool Parameter::hasSolid()
{
    return solid.number == 0 ? false : true;
}

bool Parameter::hasVirtual()
{
    return virt.number == 0 ? false : true;
}

SimTime *Parameter::getTime()
{
    return &time;
}