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

    bool print_pressure;
    int  pressure_samples;
    int  representative_cells;
    double pressure_sum;

    bool print_velocity_profile;
    int  velocity_profile_samples;

    int sample_every_n;

    CIniFile *ini;
	
    FILE *velocity_file;
    FILE *temperature_file;
    FILE *pressure_file;

    StatisticsSampler(System *system, CIniFile *ini);
	void sample();
	void calculateTemperature();
    void calculatePressure();
	void calculateVelocityProfile();
	void finish();
    double getTemperature();
    double getPressure();
};

#endif
