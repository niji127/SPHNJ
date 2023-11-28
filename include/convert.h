#ifndef _CONVERT_H
#define _CONVERT_H

#include <string>
std::string convertUnit(int size, int digit)
{
    if (size < 0)
        return "error size\n";

    int KB_size = 1024;
    int MB_size = 1024 * KB_size;
    int GB_size = 1024 * MB_size;

    if (size > GB_size)
    {
        size = pow(10.0f, (float)digit) * (float)size / (float)GB_size;
        float float_size = pow(10.0f, -(float)digit) * (float)size;
        std::string str = std::to_string(float_size);
        return str.substr(0, str.find(".") + digit + 1) + "GB";
    }
    else if (size > MB_size)
    {
        size = pow(10.0f, (float)digit) * (float)size / (float)MB_size;
        float float_size = pow(10.0f, -(float)digit) * (float)size;
        std::string str = std::to_string(float_size);
        return str.substr(0, str.find(".") + digit + 1) + "MB";
    }
    else if (size > KB_size)
    {
        size = pow(10.0f, (float)digit) * (float)size / (float)KB_size;
        float float_size = pow(10.0f, -(float)digit) * (float)size;
        std::string str = std::to_string(float_size);
        return str.substr(0, str.find(".") + digit + 1) + "KB";
    }
    else
        return std::to_string(size) + "B";
}

#endif