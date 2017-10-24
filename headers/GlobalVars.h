#ifndef _GLOBAL_VARS_H_
    #define _GLOBAL_VARS_H_

// In this header we define the global variables we will use

    // the constants of nature
    const double \
        // mass of the sun in kg
        GV_MSUN = 1.9891e30, \
        // gravitational constant in  m^3*kg^-1*s^-2
        GV_G = 6.67384e-11, \
        // speed of light in m*s^-1
        GV_C = 299792458, \
        // FIXME: What is GV_MBAR constant?
        GV_MBAR = 1.66e-24, \
        // Nice pi number
        GV_PI = 3.14159265358979323846;

    double \
    // derived units for factors, calculated at calculate_gv_units in main.c
        // density is in g/cm^3
        GV_DENSITY_UNITS, \
        // pressure is in dyns/cm^3
        GV_PRESSURE_UNITS, \
        // units for distance in km
        GV_RADCORD_UNITS, \
    // derived units for factors, calculated at calculate_gv_units in main.c
        // the machine epsilon itself
        GV_MACHINE_EPSILON, \
        // the square root of the machine epsilon
        GV_MACHINE_EPSILON_SQRT, \
        // the cube root of the machine epsilon
        GV_MACHINE_EPSILON_CUBROOT;

#endif
