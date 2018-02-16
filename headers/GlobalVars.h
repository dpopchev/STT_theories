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
    extern void dvector_copy_to_index( double *src, double *dst, int size, int *indexes);
    extern void dvector_copy_to_index_opps( double *src, double *dst, int size, int *indexes);

    extern void ode_logistics_init( ODEsystemStruct **arg );
    extern void ode_logistics_free( ODEsystemStruct **arg );
    extern void ode_logistics_compute_parameters( ODEsystemStruct *arg );

    extern FILE* open_file_to_WRITE_ResultFile( ODEsystemStruct *arg );
    extern FILE* open_file_to_APPEND_ResultFile( ODEsystemStruct *arg );
    extern FILE* open_file_to_WRITE_LivePlot( ODEsystemStruct *arg );
    extern FILE* open_file_to_APPEND_LivePlot( ODEsystemStruct *arg );
    extern FILE* open_file_to_WRITE_LivePlot_solver( ODEsystemStruct *arg );
    extern FILE* open_file_to_APPEND_LivePlot_solver( ODEsystemStruct *arg );

    // odeint vars and touches
    // they are included inside ODEsystemStruct
    extern int kmax, kount;
    extern double *xp, **yp, dxsav, *rhop, rho_tmp;
    extern void odeint(
        double ystart[], int nvar, double x1, double x2, int *nok, int *nbad,
        void (*derivs)(double, double [], double [])
    );

    extern void odeint_info_print_stdout(void);
    extern void odeint_info_print_ResultFile(FILE *fp);

    extern void newt(
        double x[], int n, int *check, void (*vecfunc)(int, double [], double [])
    );
    extern void newt_info_print_stdout(void);
    extern void newt_info_print_ResultFile(FILE *fp);

    extern void shooting_init(ShootingVarsStruct **arg);
    extern void shooting_free(ShootingVarsStruct **arg);
    extern void shooting_info_print_stdout(ShootingVarsStruct *arg);
    extern void shooting_info_print_ResultFile(ShootingVarsStruct *arg, FILE *fp);
    extern void shooting_check(ShootingVarsStruct *arg_shoot, ODEsystemStruct *arg_ode);

    extern void ode_phiScal_init( ODEsystemStruct **arg );
    extern void ode_phiScal_free( ODEsystemStruct **arg );
    extern void ode_phiScal_compute_parameters( ODEsystemStruct *arg );

    extern void eos_init(EOSmodelInfoStruct **eos);
    extern void eos_print_info_stdout(EOSmodelInfoStruct *eos);
    extern void eos_print_info_ResultFile( EOSmodelInfoStruct *eos, FILE *fp );
    extern void eos_free(EOSmodelInfoStruct **eos);
    extern void EOSeq( EOSmodelInfoStruct *eos, double pressure );

    extern double pow10(double);
#endif
