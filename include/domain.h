#ifndef _DOMAIN_H
#define _DOMAIN_H

#include "config.h"
#include "kernel.h"

class Domain
{
public:
    float cell_xmin[2], cell_xmax[2];
    float domain_xmin[2], domain_xmax[2];
    int cell_number[2];
    float interval[2];
    float interval_inv[2];
    int cell_number_total;

    bool readParameter(const Config &config);
    void initiate();
    void initiate(const Kernel &Kernel);
    bool readNumber(const std::string &path);
};

#endif