#pragma once
#include <math.h>

static double m0 = 6.63352065e-26;  // SI
static double L0 = 1e-6;            // SI
static double E0 = 1.65088e-21;     // SI
static double kb = 1.3806503e-23;   // SI

static double t0 = L0*sqrt(m0/E0);
static double F0 = E0/L0;
static double T0 = E0/kb;
static double P0 = F0/L0/L0;
static double v0 = L0/t0;
static double visc0 = P0*t0;

class UnitConverter
{
public:
    UnitConverter();

    double pressure_to_SI(double P);
    double pressure_from_SI(double P);

    double temperature_to_SI(double T);
    double temperature_from_SI(double T);

    double mass_to_SI(double m);
    double mass_from_SI(double m);

    double length_to_SI(double L);
    double length_from_SI(double L);

    double force_to_SI(double F);
    double force_from_SI(double F);

    double energy_to_SI(double E);
    double energy_from_SI(double E);

    double time_to_SI(double t);
    double time_from_SI(double t);

    double velocity_to_SI(double v);
    double velocity_from_SI(double v);

    double viscosity_to_SI(double v);
    double viscosity_from_SI(double v);
};
