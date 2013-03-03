#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

#include <stdio.h>
#include <system.h>
#include <fstream>
#include <CIniFile.h>
#include <settings.h>

class StatisticsSampler {

private:
	System *system;
    Settings *settings;
public:
    StatisticsSampler(System *system);
    void sample();
    double kinetic_energy, temperature;


};

#endif
