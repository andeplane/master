#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include <system.h>
#include <statisticssampler.h>
#include <defines.h>
#include <CIniFile.h>
#include <unitconverter.h>

using namespace std;

int main(int args, char* argv[]) {
    CIniFile ini;

    ini.load("../dsmc.ini");
    bool print_positions = ini.getbool("print_positions");
    int timesteps = ini.getint("timesteps");
    int print_every_n_step = ini.getint("print_every_n_step");

    System system;
    system.initialize(ini);
    StatisticsSampler sampler(&system, &ini);
    // StatisticsSampler *sampler = new StatisticsSampler(&system, timesteps);

    FILE *positions = 0;
    if(print_positions) {
        positions = fopen("pos.xyz","w");
        system.printPositionsToFile(positions);
    }

    for(int i=0;i<timesteps;i++) {
        if(timesteps >= 100 && !(i%(timesteps/100))) {
            printf("%d%%..",(100*i)/timesteps);
            fflush(stdout);
        }
        system.step();

        sampler.sample();

        if(print_positions && !(i%print_every_n_step)) system.printPositionsToFile(positions);
    }
    sampler.calculate_pressure();
    sampler.calculate_velocity_field();

    sampler.finish();

    UnitConverter uc;

    double dvx0 = system.dvx0;
    double dvx1 = system.dvx1;
    double P0 = system.eff_num*1.0*dvx0/system.t/system.width;
    double P1 = system.eff_num*1.0*dvx1/system.t/system.width;
    cout << dvx0 << endl;
    cout << dvx1 << endl;
    cout << P0 << endl;
    cout << P1 << endl;

    double vgrad = 2*VWALL/system.height;

    double visc = 0.5*(P1-P0)/vgrad;  // Average viscosity
    cout << "Viscosity my units: " << visc << endl;
    cout << "Viscosity factor: " << uc.viscosity_to_SI(1.0) << endl;
    cout << "Viscosity: " << uc.viscosity_to_SI(visc) << endl;
    double eta = 5.*M_PI/32.*uc.mass_to_SI(1.0)*system.density*(2./sqrt(M_PI)*uc.velocity_to_SI(system.mpv))*uc.length_to_SI(system.mfp);
    double eta_my_units = 5.*M_PI/32.*1.0*system.density*(2./sqrt(M_PI)*system.mpv)*system.mfp;
    cout << "Theoretical viscosity: " << eta << endl;
    cout << "Theoretical viscosity (my units): " << eta_my_units << endl;
    cout << "Theoretical viscosity (converted): " << uc.viscosity_to_SI(eta_my_units) << endl;

    printf("100%%\n\n");
    printf("Time consumption: \n");
    printf("Sorting: %.2f s\n",system.time_consumption[SORT]);
    printf("Moving: %.2f s\n",system.time_consumption[MOVE]);
    printf("Collisions: %.2f s\n",system.time_consumption[COLLIDE]);
    printf("Sampling: %.2f s\n",system.time_consumption[SAMPLE]);
    printf("Summary:\n");
    printf("Collisions: %d\n",system.collisions);

    printf("System volume: %f\n",system.volume);
    printf("Average temperature: %.3f\n",sampler.get_temperature());
    printf("Average pressure: %.3f\n",sampler.get_pressure());
    printf("V*P: %.3f\n",sampler.get_pressure()*system.volume);

    return 0;
}
