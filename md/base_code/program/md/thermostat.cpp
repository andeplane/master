#include <thermostat.h>
#include <system.h>
#include <statisticssampler.h>
#include <system.h>
#include <math.h>
#include <iostream>
#include <unitconverter.h>
#include <mdtimer.h>

using namespace std;

Thermostat::Thermostat(double relaxation_time_)
{
    relaxation_time = relaxation_time_;
}

void Thermostat::apply(StatisticsSampler *sampler, System *system, const double &temperature) {
    system->mdtimer->start_thermostat();
    sampler->sample_temperature();

    double berendsen_factor = sqrt(1 + system->dt/relaxation_time*(system->unit_converter->temperature_from_SI(temperature)/sampler->temperature - 1));

    for(unsigned long i=0;i<system->num_atoms_local;i++) {
        for(short a=0;a<3;a++) system->velocities[3*i+a] *= berendsen_factor;
    }
    system->mdtimer->end_thermostat();
}
