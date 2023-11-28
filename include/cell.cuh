#ifndef _CELL_H
#define _CELL_H

#include "parameter.h"
#include "particle.cuh"
#include "global.h"

#include "thrust/device_vector.h"
#include "thrust/sort.h"

class Cell
{
public:
    int neighbor_cell_id[9];

    int fluid_num;
    int solid_num;
    int virt_num;
    int particle_list[CELL_SIZE];

    int fluid_num_previous;
    int solid_num_previous;
    int virt_num_previous;
    int particle_list_previous[CELL_SIZE];

    size_t getMemoryUsed()
    {
        return sizeof(Cell);
    }

    Cell()
    {
        for (int i = 0; i < 9; i++)
            neighbor_cell_id[i] = -1;
        fluid_num = 0;
        solid_num = 0;
        virt_num = 0;
        for (int i = 0; i < CELL_SIZE; i++)
            particle_list[i] = -1;

        fluid_num_previous = 0;
        solid_num_previous = 0;
        virt_num_previous = 0;
        for (int i = 0; i < CELL_SIZE; i++)
            particle_list_previous[i] = -1;
    }
};

#endif
