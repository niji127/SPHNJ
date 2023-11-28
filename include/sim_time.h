#ifndef _SIM_TIME_H
#define _SIM_TIME_H

#include "config.h"

class SimTime
{
public:
    float dt, current_time;
    int istart, iend, i;
    int file_step;
    int solid_sub_step, shift_per_step;
    int result_interval, sort_every_step;

    bool readParameter(const Config &config);
    void initiate();
    bool addStep(); // return true if continue
    bool isFirstStep();
    bool isOutputData();
    bool isSort();
    bool readNumber(const std::string &path);
};

#endif