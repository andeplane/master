#pragma once
#include <settings.h>
#include <system.h>

class DSMC_IO
{
public:
    Settings *settings;
    System *system;

    DSMC_IO(Settings *settings_, System *system_);
    void state_to_file_binary();
    void state_from_file_binary();
    FILE *velocity_file;
    FILE *temperature_file;
    FILE *velocity_field_file_x;
    FILE *velocity_field_file_y;
    FILE *state_file;
};
