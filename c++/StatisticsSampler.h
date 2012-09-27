#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

#include <stdio.h>
#include "System.h"
#include <fstream>

class StatisticsSampler {
private:
	System *system;
public:
	bool printTemperature;
	int *numberOfSamples;

	// Viscosity sampling variables
	vec delta_v_tot;
	vec delta_v_err;
	FILE *temperatureFile;
	StatisticsSampler(System *system, int timesteps);
	void sample();
	void calculateTemperature();
	void calculateViscosity();
	void printViscosity();
};

#endif