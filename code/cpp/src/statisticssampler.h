#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

#include <stdio.h>
#include <system.h>
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

    bool print_velocity_field;

    int sample_every_n;

    CIniFile *ini;
	
    FILE *velocity_file;
    FILE *temperature_file;
    FILE *velocity_field_file_x;
    FILE *velocity_field_file_y;


    StatisticsSampler(System *system, CIniFile *ini);
	void sample();
    void calculate_temperature();
    void calculate_velocity_profile();
    void calculate_velocity_field();
    void update_cell_statistics();
    double calculate_diffusion_constant();
    mat calculate_pressure_tensor(vector<Molecule*> molecules);
    vector<mat> calculate_local_pressure_tensor();
    mat calculate_global_pressure_tensor();

	void finish();
    double get_temperature();
};

#endif
