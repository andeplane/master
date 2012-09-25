#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

#include <stdio.h>
#include "System.h"
#include <fstream>

class StatisticsSampler {
private:
	System *system;
public:
	bool temperature;
	bool pressure;
	bool energy;
	bool printVelocities;
	bool diffusionConstant;

	FILE *velocityFile;
	FILE *temperatureFile;
	FILE *pressureFile;
	FILE *energyFile;
	FILE *diffusionFile;

	StatisticsSampler(System *system);
	void sample(double t);
	void calculateTemperature(double t);
	void calculatePressure(double t);
	void calculateEnergy(double t);
	void calculateVelocities(double t);
	void calculateDiffusionConstant(double t);
};

#endif