#pragma once

class System;
class Settings;
#include <stdio.h>
#include <fstream>
#include <cinifile.h>

class StatisticsSampler {

private:
	System *system;
    Settings *settings;
    unsigned long temperature_sampled_at;
    unsigned long kinetic_energy_sampled_at;
    unsigned long velocity_distribution_sampled_at;

public:
    StatisticsSampler(System *system);
    void sample();
    void sample_kinetic_energy();
    void sample_temperature();
    void sample_velocity_distribution();
    double kinetic_energy, temperature;
};
