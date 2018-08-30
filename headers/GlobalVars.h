#ifndef _GLOBAL_VARS_H_
    #define _GLOBAL_VARS_H_

// In this header we declare the global variables we will use,
// because of the way "C" works,
// they are defined in main.c.

    extern double \
        GV_MSUN,    GV_G,   GV_C,   GV_MBAR,    GV_PI, \
        GV_DENSITY_UNITS,   GV_PRESSURE_UNITS,  GV_RADCORD_UNITS, \
        GV_MACHINE_EPSILON, GV_MACHINE_EPSILON_SQRT, \
        GV_MACHINE_EPSILON_CUBROOT, GV_J_UNITS;
    extern int GV_fadj_switch;

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
    extern int kmax, kount, did_it_go_boom;
    extern double *xp, **yp, dxsav, *rhop, rho_tmp, where_it_went_boom;
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

    extern void shooting_phiScal_init(ShootingVarsStruct **arg);
    extern void shooting_phiScal_free(ShootingVarsStruct **arg);
    extern void shooting_phiScal_info_print_stdout(ShootingVarsStruct *arg);
    extern void shooting_phiScal_info_print_ResultFile(ShootingVarsStruct *arg, FILE *fp);
    extern void shooting_phiScal_check(ShootingVarsStruct *arg_shoot, ODEsystemStruct *arg_ode);

    extern void shooting_phiScal_I_init(ShootingVarsStruct **arg);
    extern void shooting_phiScal_I_free(ShootingVarsStruct **arg);
    extern void shooting_phiScal_I_info_print_stdout(ShootingVarsStruct *arg);
    extern void shooting_phiScal_I_info_print_ResultFile(ShootingVarsStruct *arg, FILE *fp);
    extern void shooting_phiScal_I_check(ShootingVarsStruct *arg_shoot, ODEsystemStruct *arg_ode);

    extern void ode_phiScal_init( ODEsystemStruct **arg );
    extern void ode_phiScal_free( ODEsystemStruct **arg );
    extern void ode_phiScal_compute_parameters( ODEsystemStruct *arg );

    extern void ode_phiScal_I_init( ODEsystemStruct **arg );
    extern void ode_phiScal_I_free( ODEsystemStruct **arg );
    extern void ode_phiScal_I_compute_parameters( ODEsystemStruct *arg );

    extern void eos_init(EOSmodelInfoStruct **eos);
    extern void eos_print_info_stdout(EOSmodelInfoStruct *eos);
    extern void eos_print_info_ResultFile( EOSmodelInfoStruct *eos, FILE *fp );
    extern void eos_free(EOSmodelInfoStruct **eos);
    extern void EOSeq( EOSmodelInfoStruct *eos, double pressure );

    extern double pow10(double);

    extern void single_integrate_phiScal(void);
    extern void LivePlot_open(const char *ode_name, char *eos_name, double *pars);
    extern void LivePlot_append(const char *ode_name, char *eos_name, double *pars, double R,  int n);
    extern void ResultFile_open(const char *ode_name, char *eos_name, double *pars);
    extern void ResultFile_append(
      const char *ode_name, char *eos_name, double *pars,
      double p_c, double phiScal_c, double M, double AR, double rho_c,
      double delta_phiScal
    );
    extern void single_shoot_regular(void);
    extern void iterate_pressure_shoot_regular(void);
    extern void single_shoot_regular_iterate_inf(void);
    extern void iterate_pressure_inf_shoot_regular(void);
    extern void get_phiScal_cVal_infVal_regShoot(
      double *parm_vals, double pressure, EOSmodelInfoStruct *_eos,
      double *phiScal, double *inf, int _minimal_p_power
    );
    extern void get_phiScal_cVal_infVal_fitShoot(
      double *parm_vals, double pressure, EOSmodelInfoStruct *_eos,
      double *phiScal, double *inf, int _minimal_p_power
    );

    extern void single_shoot_regular_phiScal_J(void);
    extern void LivePlot_phiScal_J_open(const char *ode_name, char *eos_name, double *pars);
    extern void LivePlot_phiScal_J_append(const char *ode_name, char *eos_name, double *pars, double R,  int n);
    extern void single_shoot_regular_phiScal_J_iterate_inf(void);
    extern void single_shoot_regular_phiScal_J_iterate_inf_iterpres(void);
    extern void ResultFile_phiScal_J_open(const char *ode_name, char *eos_name, double *pars);
    extern void ResultFile_phiScal_J_append(
        const char *ode_name, char *eos_name, double *pars,
        double p_c, double phiScal_c, double M, double AR, double rho_c, double J
    );


#endif
