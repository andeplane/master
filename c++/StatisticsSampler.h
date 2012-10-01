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
	int  temperatureSamples;
	double temperatureSum;

	FILE *temperatureFile;
	StatisticsSampler(System *system, int timesteps);
	void sample();
	void calculateTemperature();
};

#endif