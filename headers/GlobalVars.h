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

    extern void dvector_copy( double *from, double *to, int size );

    extern void ode_logistics_foo( double x, double *y, double *dydx );
    extern void ode_logistics_init( ODEsystemStruct **arg );
    extern void ode_logistics_free( ODEsystemStruct **arg );
    extern void ode_logistics_integrate( ODEsystemStruct **arg );

    // odeint vars and touches
    // they are included inside ODEsystemStruct
    extern int ODEINT_SCALING_METHOD, RKQS_STEP_METHOD, kmax, kount;
    extern double *xp, **yp, dxsav;
    extern const char *ODEINT_SCALING_METHOD_DESCRIPTION[];
    extern const char *RKQS_STEP_METHOD_DESCRIPTION[];
    extern void odeint(\
        double ystart[], int nvar, double x1, double x2, double eps, double h1, \
        double hmin, int *nok, int *nbad, \
        void (*derivs)(double, double [], double []) \
    );

#endif
