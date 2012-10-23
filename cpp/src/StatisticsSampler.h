#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

#include <stdio.h>
#include <System.h>
#include <fstream>
#include <CIniFile.h>
class StatisticsSampler {
private:
	System *system;
public:
    bool print_temperature;
    int  temperature_samples;
    double temperature_sum;

    bool print_velocity_profile;
    int  velocity_profile_samples;
	
    FILE *velocity_file;
    FILE *temperature_file;

    StatisticsSampler(System *system, CIniFile &ini);
	void sample();
	void calculateTemperature();
	void calculateVelocityProfile();
	void finish();
};

#endif
