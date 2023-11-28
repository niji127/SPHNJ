#include "domain.h"

bool Domain::readNumber(const std::string &path)
{
    return false;
}

bool Domain::readParameter(const Config &config)
{
    cell_xmin[0] = config.Read("cell_xmin", 0.0f);
    cell_xmin[1] = config.Read("cell_ymin", 0.0f);
    cell_xmax[0] = config.Read("cell_xmax", 0.0f);
    cell_xmax[1] = config.Read("cell_ymax", 0.0f);

    domain_xmin[0] = config.Read("domain_xmin", 0.0f);
    domain_xmin[1] = config.Read("domain_ymin", 0.0f);
    domain_xmax[0] = config.Read("domain_xmax", 0.0f);
    domain_xmax[1] = config.Read("domain_ymax", 0.0f);

    if (cell_xmax[0] - cell_xmin[0] == 0.0f || cell_xmax[1] - cell_xmin[1] == 0.0f)
    {
        std::cout << "failed to read cell information" << std::endl;
        return true;
    }
    if (domain_xmax[0] - domain_xmin[0] == 0.0f || domain_xmax[1] - domain_xmin[1] == 0.0f)
    {
        std::cout << "failed to read domain information" << std::endl;
        return true;
    }

    return false;
}

void Domain::initiate()
{}

void Domain::initiate(const Kernel &kernel)
{
    float smoothing_length = kernel.smoothing_length;
    float cell_size_max = kernel.impact_length_hsml_ratio * smoothing_length;

    for (int i = 0; i < 2; i++)
    {
        cell_number[i] = (int)((cell_xmax[i] - cell_xmin[i]) / cell_size_max);
        interval[i] = (cell_xmax[i] - cell_xmin[i]) / (float)cell_number[i];
        interval_inv[i] = 1.0f / interval[i];
    }
    cell_number_total = cell_number[0] * cell_number[1];
}