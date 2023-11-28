#include "physics.h"

void Physics::initiate() {}
bool Physics::readNumber(const std::string &path)
{
    return false;
}

bool Physics::readParameter(const Config &config)
{
    {
        float gravity_y = config.Read("gravity", 9.81f);
        gravity[0] = 0.0f;
        gravity[1] = -gravity_y;

        bool transient_type_bool = false;
        transient_type_bool = config.Read("is_steady", false);
        if (transient_type_bool)
            trasient_type = STEADY;
        else
            trasient_type = TRANSIENT;

        return false;
    }
}