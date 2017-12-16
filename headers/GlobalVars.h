#ifndef _GLOBAL_VARS_H_
    #define _GLOBAL_VARS_H_

// In this header we declare the global variables we will use,
// because of the way "C" works,
// they are defined in main.c.

    extern double \
        GV_MSUN,    GV_G,   GV_C,   GV_MBAR,    GV_PI, \
        GV_DENSITY_UNITS,   GV_PRESSURE_UNITS,  GV_RADCORD_UNITS, \
        GV_MACHINE_EPSILON, GV_MACHINE_EPSILON_SQRT, \
        GV_MACHINE_EPSILON_CUBROOT;
//    extern char *ODEsystem_names[];
//    extern char *ODEsystem_vars[];

#endif
