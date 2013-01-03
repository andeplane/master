#pragma once

class UnitConverter
{
private:
    void calculate_units();
public:
    UnitConverter();

    static double m0 = 6.63352065e-26;  // SI
    static double L0 = 1e-6;            // SI
    // static double E0 = 1.0318e-2;    // eV
    static double E0 = 1.65088e-21;     // SI
    static double kb = 1.3806503e-23;   // SI

    // Will be set in calculate_units()
    static double t0 = 1.0;
    static double F0 = 1.0;
    static double T0 = 1.0;
    static double PA0 = 1.0;
};
