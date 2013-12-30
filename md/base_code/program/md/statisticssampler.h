#pragma once

#include <stdio.h>
#include <system.h>
#include <fstream>

class System;
class Settings;

class StatisticsSampler {
private:
    System *system;
    Settings *settings;
    unsigned long pressure_sampled_at;
    unsigned long temperature_sampled_at;
    unsigned long kinetic_energy_sampled_at;
    unsigned long potential_energy_sampled_at;

public:
    StatisticsSampler(System *system);
    void sample();
    void sample_temperature();
    void sample_kinetic_energy();
    void sample_potential_energy();
    void sample_pressure();
    void sample_velocity_distribution();
    void sample_momentum_cm();
    double pressure;
    double kinetic_energy;
    double potential_energy;
    double temperature;
    double v_cm[3];
    double p_cm[3];
};
