#include "ExternalHeaders.h"

#define DEBUGGING_ode_phiScal_I_foo 0
#define DEBUGGING_ode_phiScal_I_init 0
#define DEBUGGING_ode_phiScal_I_free 0
#define DEBUGGING_ode_phiScal_I_integrate 0
#define DEBUGGING_ode_phiScal_I_compute_parameters 0
#define DEBUGGING_ode_phiScal_I_change_central_value 0

// General variables
#define ODE_VARS_NAME_LENGHT 16
#define MAX_ARR_SIZE 16

// ODE system general info
#define ODE_INDEX 2
#define ODE_DESCRIPTION "STT_phiScal_I"

// ODE count of equations and initial variable per equation
#define ODE_EQS_COUNT 9
#define ODE_Y_INIT_VAL 0, 0, 5e-4, 0, 0, 0, 0 ,0 ,0

// ODE names/symbols for independent and for each dependent variable
#define ODE_NAME_INDEP "r"
#define ODE_NAME_DEP \
  "phiScal", "Q", "p", "LambdaMetr", "m", "PhiMetr", "Z", "Omega", "J"

// odeint mount of points we want to print if any
#define ODE_POINTS_COUNT 1e4

// amount of ode free parameters, name/symbol for each one
#define ODE_FREE_PARM_COUNT_ALL 3
#define ODE_FREE_PARM_NAME "beta", "m", "lambda"

// beta
// beta to be up to -10
// 3
#define ODE_FREE_PARM_COUNT_1 1
// 0, -6, -10
#define ODE_FREE_PARM_VALS_1 -6

// m
// 3
#define ODE_FREE_PARM_COUNT_2 1
// 0, 1e-3, 5e-2
#define ODE_FREE_PARM_VALS_2 5e-2

// lambda
// 5
#define ODE_FREE_PARM_COUNT_3 1
// 0, 1e-1, 1, 10, 100
#define ODE_FREE_PARM_VALS_3 0

// the interval for the independent variable
#define ODE_INDEP_INIT 1e-30
#define ODE_INDEP_FINAL 1e2

// if we are to change vary any of the initial values of the ODE
// index, start and end value, also the step to achieve it
#define Y_INDEX_CHANGE 3
#define INITIAL_Y_START 2.65e-4 // 1e-5
#define INITIAL_Y_END 5e-3 // 5e-3
#define INITIAL_Y_STEP 1 // 1e-5

// if we have greater difference than the power of this we stop iterating over guesses
#define MIN_DIFF_POW -3
// if the guessed value for scalar field is less than -8 the re is no point in trying
#define MIN_PHISCAL_POW -3
// the coefficient to increase the infinity
#define INCREASE_COEF 0.3

// how many 1/m we will try the fitting
#define FITTING_MASSIVE_COEF 1

static double *GV_PARAMETERS_VALUES, AR, R, tiny = 1e-30, rho_c;
static EOSmodelInfoStruct *eos;
static int too_small_to_count_pressure_pow = -16, pressure_index = 3;

static void ode_phiScal_I_foo(double x, double *y, double *dydx){

    //for(int i=1; i <= 9; i++){
        //if(!isfinite(dydx[i])){
            //printf("\n at %e y[ %d ] not finite ---> %e \n", x,i, y[i]);
        //}
    //}

    double \
        beta_phiScal = GV_PARAMETERS_VALUES[1],
        m_phiScal = GV_PARAMETERS_VALUES[2],
        lambda_phiScal = GV_PARAMETERS_VALUES[3],
        r = x,
        phiScal = y[1],
        Q = y[2],
        p = y[3],
        LambdaMetr = y[4],
        m = y[5],
        PhiMetr = y[6],
        Z = y[7],
        Omega = y[8],
        J = y[9],
        \
        Vhat = 2*pow(m_phiScal,2)*pow(phiScal,2) + lambda_phiScal*pow(phiScal,4),
        Vhat_dphiScal = 4*pow(m_phiScal,2)*phiScal + 4*lambda_phiScal*pow(phiScal,3),
        alpha = beta_phiScal*phiScal,
        A = exp(1.0/2 * beta_phiScal * pow(phiScal,2)),
        exp_2LambdaMetr = exp(2*LambdaMetr ),
        \
        step4_A = pow(A,4),
        rho;

    if( !R && p && floor(log10(fabs( p + tiny))) <= too_small_to_count_pressure_pow ){
        R = r;
        AR = A*r;
        p = 0;
        rho = 0;
    }else if( floor(log10(fabs( p + tiny))) <= too_small_to_count_pressure_pow ){
        p = 0;
        rho = 0;
    }else{
        EOSeq(eos,p);
        rho = eos->current;

        if(!rho_c){
            rho_c = rho;
        }
    }

    rho_tmp = rho;

    double \
      Omega_dr = Z, \
      \
      PhiMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*p*exp_2LambdaMetr
          + pow(Q,2)
          - 1.0/2*exp_2LambdaMetr*Vhat
          + pow(1/r,2)*(exp_2LambdaMetr-1)
        ), \
      \
      phiScal_dr = Q, \
      \
      LambdaMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*rho*exp_2LambdaMetr
          + pow(Q,2)
          + 1.0/2*exp_2LambdaMetr*Vhat
          - pow(1/r,2)*(exp_2LambdaMetr-1)
        ), \
      \
      Q_dr = \
        4*GV_PI*alpha*step4_A*(rho - 3*p)*exp_2LambdaMetr \
        + 1.0/4*Vhat_dphiScal*exp_2LambdaMetr \
        - (PhiMetr_dr - LambdaMetr_dr + 2/r )*Q, \
      \
      p_dr = -( rho + p ) * ( PhiMetr_dr + alpha*Q ), \
      \
      m_dr = \
        pow(r,2) * (
          4 * GV_PI * step4_A* rho
          + 1.0/2 * 1/exp_2LambdaMetr * pow(Q,2)
          + 1.0/4 * Vhat
        ), \
      \
      Z_dr = \
        16*GV_PI*step4_A*(rho+p)*exp_2LambdaMetr*Omega \
        - Z * ( 4/r - (PhiMetr_dr + LambdaMetr_dr) ), \
      \
      J_dr = 8.0/3*GV_PI*step4_A*(rho+p)*exp(LambdaMetr-PhiMetr)*pow(r,4)*Omega;

    dydx[1] = phiScal_dr;
    dydx[2] = Q_dr;
    dydx[3] = p_dr;
    dydx[4] = LambdaMetr_dr;
    dydx[5] = m_dr;
    dydx[6] = PhiMetr_dr;
    dydx[7] = Z_dr;
    dydx[8] = Omega_dr;
    dydx[9] = J_dr;

    //for(int i=1; i <= 9; i++){
        //if(!isfinite(dydx[i])){
            //printf("\n at %e dydx[ %d ] not finite ---> %e \n", x,i, dydx[i]);
        //}
    //}

    return;
}

void ode_phiScal_I_init( ODEsystemStruct **arg ){

    const char \
        function_path[] = "ode_phiScal_I.c ode_phiScal_I_init: ", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_I_init){
        printf("%s %s starting \n", identation, function_path);
    }

    if( ( *arg = calloc(1,sizeof(ODEsystemStruct)) ) == NULL ){
        printf("%s %s error in calloc for ode, line 92 \n", identation, function_path);
        exit(92);
    };

    // ODE system index
    (*arg)->index = ODE_INDEX;

    if(DEBUGGING_ode_phiScal_I_init){
        printf( \
          "%s %s ODE system index %d \n", \
          identation, function_path, (*arg)->index \
        );
    };

    // total number of ODE equations
    (*arg)->eqs_count = ODE_EQS_COUNT;
    if(!(*arg)->eqs_count){
        printf(\
          "%s %s HOW COULD IT BE ODE SYSTEM WITH SIZE 0!!!!, line 101", \
          identation,function_path \
        );

        exit(101);
    }else if(DEBUGGING_ode_phiScal_I_init){
        printf( \
          "%s %s ODE system eqs count %d \n", \
          identation, function_path, (*arg)->eqs_count \
        );
    };

    // initial values for *y
    // they should be as many as (*arg)->eqs_count
    // the first value is not counted, since dvector starts from 1 !!!
    double y[] = {
        1e10, // THIS VALUE DOES NOT TAKE ANY EFFECT, since dvector starts from index 1
        ODE_Y_INIT_VAL // initial population
    };

    (*arg)->y = dvector(1, (*arg)->eqs_count);
    dvector_copy(y, (*arg)->y, (*arg)->eqs_count);

    if(DEBUGGING_ode_phiScal_I_init){
        for(int i=1; i <= (*arg)->eqs_count; i++){
            printf(\
              "%s %s y[%d] = %.3e \n", \
              identation,function_path, i, (*arg)->y[i] \
            );
        }
    }

    // index for the y whos initial value will be changed
    // if any, if not just 0
    (*arg)->index_of_y_to_change = Y_INDEX_CHANGE;

    if((*arg)->index_of_y_to_change){
        (*arg)->initial_y_start = INITIAL_Y_START;
        (*arg)->initial_y_end = INITIAL_Y_END;
        (*arg)->initial_y_current = (*arg)->initial_y_start;
        (*arg)->initial_y_step = INITIAL_Y_STEP;

        if(DEBUGGING_ode_phiScal_I_init){

            printf(
              "%s %s We will change initial value of y[%d] from %.3e to %.3e with step %.3e\n",
              identation,function_path,
              (*arg)->index_of_y_to_change,
              (*arg)->initial_y_start, (*arg)->initial_y_end, (*arg)->initial_y_step
            );
        }
    }else{
        (*arg)->initial_y_start = 0;
        (*arg)->initial_y_end = 0;
        (*arg)->initial_y_current = 0;

        if(DEBUGGING_ode_phiScal_I_init){
            printf("%s %s We will NOT change \n",identation,function_path);
        }
    }

    // variable names per index
    // the 0 index is the name of the independent variable
    const char *ode_names[] = {
        ODE_NAME_INDEP, // independent
        ODE_NAME_DEP // dependent
    };

    (*arg)->name_vars = calloc( (size_t)(*arg)->eqs_count+1, sizeof(char*));

    for(int i=0; i <=(*arg)->eqs_count; i++ ){
        (*arg)->name_vars[i] = calloc( ODE_VARS_NAME_LENGHT, sizeof(char));

        strcpy((*arg)->name_vars[i],ode_names[i]);

        if(DEBUGGING_ode_phiScal_I_init){
            i ? \
                printf(\
                  "%s %s the dependent [%d] variable is %s \n", \
                  identation,function_path, i,(*arg)->name_vars[i] \
                ) \
                : \
                printf(\
                  "%s %s the independent variable is %s \n", \
                  identation,function_path, (*arg)->name_vars[i] \
                ) \
            ;
        }
    }

    // name name, index and simple description
    strcpy((*arg)->name_system, ODE_DESCRIPTION);

    if(DEBUGGING_ode_phiScal_I_init){
        printf(\
          "%s %s System name and description: %s \n", \
          identation,function_path, (*arg)->name_system \
        );
    }

    // amount of points we want to save and print, if any
    (*arg)->points_count = ODE_POINTS_COUNT;
    if( (*arg)->points_count){

        kmax = (*arg)->points_count;

        (*arg)->points_x = dvector(1, (*arg)->points_count);
        (*arg)->points_y = dmatrix(1, (*arg)->eqs_count, 1, (*arg)->points_count);
        (*arg)->points_rho = dvector(1, (*arg)->points_count);

        xp = (*arg)->points_x;
        yp = (*arg)->points_y;
        rhop = (*arg)->points_rho;

    }else{

        xp = NULL;
        yp = NULL;
        rhop = NULL;

    }

    if(DEBUGGING_ode_phiScal_I_init){
        printf(\
          "%s %s We want to print %d points \n", \
          identation,function_path, (*arg)->points_count \
        );
    }

    // count OK steps done by the integrator
    (*arg)->nok = 0;

    if(DEBUGGING_ode_phiScal_I_init){
        printf(\
          "%s %s We have %d OK points \n", \
          identation,function_path, (*arg)->nok \
        );
    }

    // count BAD steps done by the integrator
    (*arg)->nbad = 0;

    if(DEBUGGING_ode_phiScal_I_init){
        printf(\
          "%s %s We have %d BAD points \n", \
          identation,function_path, (*arg)->nok \
        );
    }

    // how many free parameters the system has
    (*arg)->free_parmeters_count_all = ODE_FREE_PARM_COUNT_ALL;
    if(DEBUGGING_ode_phiScal_I_init){
        printf(\
          "%s %s We have %d free parameters in the system \n", \
          identation,function_path, (*arg)->free_parmeters_count_all \
        );
    }

    const char *free_parameters_names[] = {
        "1e10", // empty since it will not be red
        ODE_FREE_PARM_NAME // parameter names
    };

    (*arg)->free_parmeters_names = \
        calloc( (size_t)(*arg)->free_parmeters_count_all + 1, sizeof(char*) );

    for(int i=1; i <=(*arg)->free_parmeters_count_all; i++ ){

        (*arg)->free_parmeters_names[i] = calloc( ODE_VARS_NAME_LENGHT, sizeof(char));

        strcpy((*arg)->free_parmeters_names[i],free_parameters_names[i]);

        if(DEBUGGING_ode_phiScal_I_init){
            printf(\
              "%s %s [%d] free parameter name is %s\n", \
              identation,function_path, i, (*arg)->free_parmeters_names[i]\
            );
        }
    }

    (*arg)->free_parmeters_values = dvector(1, (*arg)->free_parmeters_count_all);
    GV_PARAMETERS_VALUES = (*arg)->free_parmeters_values;

    if(DEBUGGING_ode_phiScal_I_init){
        printf(
          "%s %s The free parameters of the system are at addresses:\n \t \
          local variable %p; part of the ode system %p \n",
          identation,function_path,
          (void*)GV_PARAMETERS_VALUES, (void*)(*arg)->free_parmeters_values
        );
    }

    (*arg)->free_parameters_count_each = \
        calloc( (size_t)(*arg)->free_parmeters_count_all+1, sizeof(int));

    (*arg)->free_parmeters_values_all = \
        calloc( (size_t)(*arg)->free_parmeters_count_all+1, sizeof(double*));

    // all values for the free parameters
    double free_parmeters_values_all[][MAX_ARR_SIZE] = {
        { 1e10 },  // index starts from 1 so this will not be taken into account
        { 1e10, ODE_FREE_PARM_VALS_1 },  // beta
        { 1e10, ODE_FREE_PARM_VALS_2 },  // m
        { 1e10, ODE_FREE_PARM_VALS_3 }  // lambda
    };

    for(int i=1; i <= (*arg)->free_parmeters_count_all; i++){
        (*arg)->free_parmeters_values[i] = free_parmeters_values_all[i][1];

        if(DEBUGGING_ode_phiScal_I_init){
            printf(\
              "%s %s Initial value for free parameter [%d]->%s = %.3e \n", \
              identation,function_path, i, \
              (*arg)->free_parmeters_names[i], (*arg)->free_parmeters_values[i] \
            );
        }
    }

    // how many values of each parameter we will investigate
    int free_parameters_count_each[] = {
        1e5, // index starts from 1 so this will not be taken into account
        ODE_FREE_PARM_COUNT_1,// count all values of parameter 1 we want to check
        ODE_FREE_PARM_COUNT_2, // count all values of parameter 1 we want to check
        ODE_FREE_PARM_COUNT_3 // count all values of parameter 1 we want to check
    };

    for(int i=1; i <= (*arg)->free_parmeters_count_all; i++){

        (*arg)->free_parameters_count_each[i] = \
            free_parameters_count_each[i];

        (*arg)->free_parmeters_values_all[i] = \
            calloc( (size_t)(*arg)->free_parameters_count_each[i]+1, sizeof(double));

        for(int l=1; l <= (*arg)->free_parameters_count_each[i]; l++){

            (*arg)->free_parmeters_values_all[i][l] = \
                free_parmeters_values_all[i][l];

            if(DEBUGGING_ode_phiScal_I_init){
                printf(\
                  "%s %s Possible value %d/%d for parameter [%d]-> %s = %.3e \n", \
                  identation,function_path, l, free_parameters_count_each[i],i, \
                  (*arg)->free_parmeters_names[i], (*arg)->free_parmeters_values_all[i][l] \
                );
            }
        }
    }

    // we will integrate starting from this point
    (*arg)->x_initial = ODE_INDEP_INIT;

    // until we reach this one
    (*arg)->x_final = ODE_INDEP_FINAL;

    if(DEBUGGING_ode_phiScal_I_init){
        printf(\
          "%s %s Independent variable %s initial %.3e and final %.3e value \n", \
          identation,function_path, (*arg)->name_vars[0], \
          (*arg)->x_initial, (*arg)->x_final \
        );
    }

    (*arg)->foo = &ode_phiScal_I_foo;

    if(DEBUGGING_ode_phiScal_I_init){
        printf("%s %s ending \n", identation, function_path);
    }

    return;
}

static void ode_phiScal_I_info_print_stdout( ODEsystemStruct *arg ){

    printf(
        "\n\n System info: \n\n"
        "\t General info: \n"
        "\t\t index: %d \n"
        "\t\t name: %s \n",
        arg->index,
        arg->name_system
    );

    printf(
        "\n\t Independent variable info: \n"
        "\t\t name: %s \n"
        "\t\t initial value: %.3e \n"
        "\t\t final value: %.3e \n"
        "\t\t increase inf coef: %.3e \n"
        "\t\t min phiScal power vars: %d \n"
        "\t\t min diff power for shooting vars: %d \n"
        "\t\t xfitting massive coef: %d \n",
        arg->name_vars[0], arg->x_initial, arg->x_final,
        INCREASE_COEF, MIN_PHISCAL_POW, MIN_DIFF_POW,
        FITTING_MASSIVE_COEF
    );

    printf(
        "\n\t Dependent variables info: \n"
    );
    for(int i=1; i<=arg->eqs_count; i++){
        printf(
            "\t\t initial value for %s = %.3e \n",
            arg->name_vars[i], arg->y[i]
        );
    }

    if(arg->index_of_y_to_change){
        printf(
          "\n\t\t Initial value for %s "
          "will change from %.3e to %.3e with step %.3e \n",
          arg->name_vars[arg->index_of_y_to_change],
          arg->initial_y_start,
          arg->initial_y_end,
          arg->initial_y_step
        );
    }

    printf("\n\t Free parameters info: \n");
    for(int i=1; i<=arg->free_parmeters_count_all; i++){

        printf(
            "\t\t current value for %s = %.3e \n",
            arg->free_parmeters_names[i], arg->free_parmeters_values[i]
        );

        printf("\t\t rest values for %s = ", arg->free_parmeters_names[i] );

        for(int l=1; l<=arg->free_parameters_count_each[i]; l++){
            printf(" %.3e ", arg->free_parmeters_values_all[i][l] );
        }

        printf("\n");
    }

    if(arg->points_count){
        printf("\n\t\t Integrator will record at most %d points \n", arg->points_count);
    }

    return;
}

static void ode_phiScal_I_info_print_ResultFile( ODEsystemStruct *arg){

     FILE *fp = open_file_to_APPEND_ResultFile(arg);

    fprintf(
        fp,
        "\n\n System info: \n\n"
        "\t General info: \n"
        "\t\t index: %d \n"
        "\t\t name: %s \n",
        arg->index,
        arg->name_system
    );

    fprintf(
        fp,
        "\n\t Independent variable info: \n"
        "\t\t name: %s \n"
        "\t\t initial value: %.3e \n"
        "\t\t final value: %.3e \n"
        "\t\t increase inf coef: %.3e \n"
        "\t\t min phiScal power vars: %d \n"
        "\t\t min diff power for shooting vars: %d \n"
        "\t\t xfitting massive coef: %d \n",
        arg->name_vars[0], arg->x_initial, arg->x_final,
        INCREASE_COEF, MIN_PHISCAL_POW, MIN_DIFF_POW,
        FITTING_MASSIVE_COEF
    );

    fprintf(
        fp,
        "\n\t Dependent variables info: \n"
    );
    for(int i=1; i<=arg->eqs_count; i++){
        fprintf(
            fp,
            "\t\t initial value for %s = %.3e \n",
            arg->name_vars[i], arg->y[i]
        );
    }

    if(arg->index_of_y_to_change){
        fprintf(
          fp,
          "\n\t\t Initial value for %s "
          "will change from %.3e to %.3e with step %.3e \n",
          arg->name_vars[arg->index_of_y_to_change],
          arg->initial_y_start,
          arg->initial_y_end,
          arg->initial_y_step
        );
    }

    fprintf( fp,"\n\t Free parameters info: \n");
    for(int i=1; i<=arg->free_parmeters_count_all; i++){

        fprintf(
            fp,
            "\t\t current value for %s = %.3e \n",
            arg->free_parmeters_names[i], arg->free_parmeters_values[i]
        );

        fprintf(
            fp,"\t\t rest values for %s = ", arg->free_parmeters_names[i]
        );

        for(int l=1; l<=arg->free_parameters_count_each[i]; l++){
            fprintf(
                fp," %.3e ", arg->free_parmeters_values_all[i][l]
            );
        }

        fprintf(
            fp,
            "\n"
        );
    }

    if(arg->points_count){
        fprintf(
            fp,
            "\n\t\t Integrator will record at most %d points \n",
            arg->points_count
        );
    }

    fclose(fp);

    return;
}

static void ode_phiScal_I_ResultFile_open(ODEsystemStruct *arg){

    FILE *fp = open_file_to_WRITE_ResultFile(arg);

    fclose(fp);

    return;
}

static void ode_phiScal_I_ResultFile_append_shootregular(
  ODEsystemStruct *arg, ShootingVarsStruct *shoot_vars
){

    FILE *fp = open_file_to_APPEND_ResultFile(arg);

    fprintf(
      fp,
      "%e %e %e %e %e\n",
      arg->initial_y_current,
      shoot_vars->newt_v[1],
      arg->y[5],
      AR,
      fabs(arg->y[1] - shoot_vars->known_right_values[1])
    );

    fclose(fp);

    return;
}

static void ode_phiScal_I_ResultFile_append_shootfitting(
  ODEsystemStruct *arg, ShootingVarsStruct *shoot_vars
){

    FILE *fp = open_file_to_APPEND_ResultFile(arg);

    fprintf(
      fp,
      "%e %e %e %e %e %e %e %e %e %e %e\n",
      arg->initial_y_current,
      shoot_vars->newt_v[1],
      shoot_vars->newt_v[2],
      shoot_vars->newt_v[3],
      //arg->y[5],
      //M,
      shoot_vars->newt_v[8],
      AR,
      rho_c,
      shoot_vars->newt_v[7],
      shoot_vars->newt_v[4],
      shoot_vars->newt_v[5],
      shoot_vars->newt_v[6]
    );

    fclose(fp);

    return;
}

static void ode_phiScal_I_LivePlot_open_solver(ODEsystemStruct *arg){

    FILE *fp = open_file_to_WRITE_LivePlot_solver(arg);

    fclose(fp);

    return;
}

static void ode_phiScal_I_LivePlot_append_solver(ODEsystemStruct *arg){

    FILE *fp = open_file_to_APPEND_LivePlot_solver(arg);

    if(arg->points_count){
        fprintf(
          fp,
          "# phiScal_c = %e, R = %e, phiScal_inf = %e",
          arg->points_y[1][1], R, arg->phiScal_inf
        );

        for(int i=1; i <= arg->points_count && arg->points_x[i]; i++){
            fprintf(fp,"\n%e ", arg->points_x[i]);

            for(int j=1; j <= arg->eqs_count; j++){
                fprintf(fp,"%e ", arg->points_y[j][i]);
            }

            fprintf(fp,"%e ", arg->points_rho[i]);
        }
        fprintf(fp,"\n");
    }else{
        fprintf(fp,"# NO POINTS OOOO");
    }

    //printf("\n\t\t\t Small sleep when appending to live plot %d \n", GV_fadj_switch);
    //sleep(1);

    fclose(fp);

    //sleep(1);

    return;
}

void ode_phiScal_I_free( ODEsystemStruct **arg ){

    const char \
        function_path[] = "ODE_logistic.c ode_logistics_free", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_I_free){
        printf("%s %s starting \n", identation, function_path);
    }

    free((*arg)->y);

    for(int i=0; i <=(*arg)->eqs_count; i++ ){
        free((*arg)->name_vars[i]);
    }

    free((*arg)->name_vars);

    if( (*arg)->points_count ){

        free( (*arg)->points_x);

        free_dmatrix(
          (*arg)->points_y, 1, (*arg)->eqs_count,
          1, (*arg)->points_count
        );

    }else{
        free((*arg)->points_x);
        free((*arg)->points_y);
    }

    for(int i=1; i <=(*arg)->free_parmeters_count_all; i++ ){
        free((*arg)->free_parmeters_names[i]);
    }

    free((*arg)->free_parmeters_names);

    free(GV_PARAMETERS_VALUES);

    free((*arg)->free_parameters_count_each );

    for(int i=1; i <= (*arg)->free_parmeters_count_all; i++){
        free((*arg)->free_parmeters_values_all[i]);
    }

    free((*arg)->free_parmeters_values_all);

    free((*arg)->points_rho);

    free((*arg)->eoseq_name);

    free(*arg);

    return;
}

static void ode_phiScal_I_integrate( ODEsystemStruct *arg ){

    const char function_path[] = "ode_phiScal_I.c ode_phiScal_I_integrate", \
               identation[] = "\n";

    //if(DEBUGGING_ode_phiScal_I_integrate){
        //printf("\n\n %s %s starting \n", identation, function_path);
    //}

    if(kmax){
        for(int i=1; i <= arg->points_count && xp[i]; i++){
            xp[i] = 0;
            for(int l=1; l <= arg->eqs_count; l++){
                yp[l][i] = 0;
            }
            rhop[i] = 0;
        }
    }else{
        xp = NULL;
        yp = NULL;
        rhop = NULL;
    }

    arg->nok = arg->nbad = 0;

    if(DEBUGGING_ode_phiScal_I_integrate && !GV_fadj_switch){
        printf("\n\n");
        printf("r = %e", arg->x_initial);
        for(int i=1; i<=arg->eqs_count;i++){
            printf(" %s_init = %e ", arg->name_vars[i], arg->y[i]);
        }
    }

    if(arg->x_final > arg->x_initial){
        R =  0;
        AR = 0;
        rho_c = 0;
    }

    arg->did_it_go_boom = 0;
    arg->where_it_went_boom = 0;

    //did_it_go_boom = &arg->did_it_go_boom;
    //where_it_went_boom = &arg->where_it_went_boom;

    odeint(
      arg->y, arg->eqs_count,
      arg->x_initial, arg->x_final,
      &arg->nok, &arg->nbad,
      arg->foo
    );

    //M = arg->y[5];

    if(DEBUGGING_ode_phiScal_I_integrate && !GV_fadj_switch){
        printf("\n\n");
        printf("r = %e", arg->x_final);
        for(int i=1; i<=arg->eqs_count;i++){
            printf(" %s_final = %e ", arg->name_vars[i], arg->y[i]);
        }
    }

    //if(DEBUGGING_ode_phiScal_I_integrate){
        //printf("\n%s %s ending \n\n", identation, function_path);
    //}

    //sleep(2);

    return;
}

static int \
    guess_left_n, *guess_left_indexes, \
    guess_right_n, *guess_right_indexes, \
    shoot_free_n, *shoot_free_indexes, \
    solver_try;

static double \
  *known_left_values, *known_right_values, shoot_fiting_point;

static ODEsystemStruct *guess_to_integrate;

static void ode_phiScal_I_shooting_phiScal_I_regular(int n, double *v, double *f){

    double \
      *tmp_y = dvector(1, guess_to_integrate->eqs_count),
      *tmp_f1 = dvector(1, guess_left_n );

    if(n != guess_left_n){
        printf(
          "ode_phiScal_I.c ode_phiScal_I_shooting_phiScal_I_regular something"
          "verrryyy fissshyyy %d %d line 822",
          n, guess_left_n
        );

        exit(826);
    }

    tmp_y = dvector(1, guess_to_integrate->eqs_count);
    dvector_copy(guess_to_integrate->y, tmp_y ,guess_to_integrate->eqs_count);

    dvector_copy_to_index(
      v,
      guess_to_integrate->y,
      guess_left_n,
      guess_left_indexes
    );

    dvector_copy_to_index(
      known_left_values,
      guess_to_integrate->y,
      guess_right_n,
      guess_right_indexes
    );

    ode_phiScal_I_integrate(guess_to_integrate);

    dvector_copy_to_index_opps(
      guess_to_integrate->y,
      tmp_f1,
      guess_left_n,
      guess_left_indexes
    );

    for(int i=1; i <= n; i++){
        f[i] = known_right_values[i] - tmp_f1[i];
    }

    dvector_copy(tmp_y, guess_to_integrate->y ,guess_to_integrate->eqs_count);
    free(tmp_y);
    free(tmp_f1);

    return;
}

// Not up to date !!!!
static void ode_phiScal_I_shoot_regular_execute(
  ODEsystemStruct *arg, ShootingVarsStruct *shoot_vars
){

    guess_left_n = shoot_vars->UNknown_left_n;
    guess_left_indexes = shoot_vars->UNknown_left_indexes;
    known_right_values = shoot_vars->known_right_values;

    guess_right_n = shoot_vars->UNknown_right_n;
    guess_right_indexes = shoot_vars->UNknown_right_indexes;
    known_left_values = shoot_vars->known_left_values;

    shoot_free_n = shoot_vars->shoot_free_n;
    shoot_free_indexes = shoot_vars->shoot_free_indexes;

    shoot_vars->newt_n = shoot_vars->UNknown_left_n;

    guess_to_integrate = arg;

    dvector_copy(
      shoot_vars->UNknown_left_values,
      shoot_vars->newt_v,
      shoot_vars->newt_n
    );

    //ode_phiScal_I_LivePlot_open_solver(arg);

    int newt_checker;
    newt(
      shoot_vars->newt_v,
      shoot_vars->newt_n,
      &newt_checker,
      &ode_phiScal_I_shooting_phiScal_I_regular
    );

    dvector_copy_to_index(
      shoot_vars->newt_v,
      arg->y,
      guess_left_n,
      guess_left_indexes
    );

    ode_phiScal_I_integrate(arg);

    printf(
      "\n REGULAR try p_c = %.3e; \t phi_c = %.3e -> phi_g = %.3e \t newt_check %d",
      arg->initial_y_current,
      shoot_vars->UNknown_left_values[1],
      shoot_vars->newt_v[1],
      newt_checker
    );

    //ode_phiScal_I_ResultFile_append_shootregular(arg,shoot_vars);

    guess_left_indexes = NULL;
    guess_right_indexes = NULL;
    shoot_free_indexes = NULL;
    known_left_values = NULL;
    known_right_values = NULL;
    guess_to_integrate = NULL;

    return;
}

static void ode_phiScal_I_shooting_phiScal_I_fitting(int n, double *v, double *f){

    if(n != guess_left_n + guess_right_n){
        printf(
          "\n ode_phiScal_I.c ode_phiScal_I_shooting_phiScal_I_fitting something"
          "verrryyy fissshyyy %d %d line 853 \n",
          n, guess_left_n + guess_right_n
        );
        exit(856);
    }

    double \
        *tmp_y = dvector(1, guess_to_integrate->eqs_count),
        *tmp_f1 = dvector(1, guess_left_n + guess_right_n ),
        *tmp_f2 = dvector(1, guess_left_n + guess_right_n),
        tmp_x_left = guess_to_integrate->x_initial,
        tmp_x_right = guess_to_integrate->x_final,
        coef = 0.1, coef_step = 0.05;

    dvector_copy(guess_to_integrate->y, tmp_y ,guess_to_integrate->eqs_count);

    int iterate = 1;

    // we will loop through the process of evaluating the discrepancy at the fitting point
    // by integrating from the left end in the interval of interest to a fitting point
    // and from right end of the interval of interest to a fitting point
    // if for any circumstance we have encountered an "boom" solution
    // we change the fitting point with little bit from left/right of the
    // point of "explosion" depending on whether we are left ot fit or right to fit
    while(iterate){

        dvector_copy_to_index(
          v,
          guess_to_integrate->y,
          guess_left_n,
          guess_left_indexes
        );

        dvector_copy_to_index(
          known_left_values,
          guess_to_integrate->y,
          guess_right_n,
          guess_right_indexes
        );

        guess_to_integrate->x_initial = tmp_x_left;
        guess_to_integrate->x_final = shoot_fiting_point;
        guess_to_integrate->phiScal_inf = shoot_fiting_point;

        ode_phiScal_I_integrate(guess_to_integrate);

        // append only if we are not calculating fadj in the solver
        if(!GV_fadj_switch){
            ode_phiScal_I_LivePlot_append_solver( guess_to_integrate );
        }

        // we just make sure the pressure is not 0 for the integration right to left
        guess_to_integrate->y[pressure_index] = 0;

        if(guess_to_integrate->did_it_go_boom){

            printf(
              "\n\t\t\t\t\t it whent booommmmm left to right at %e when using fitt %e \n",
              guess_to_integrate->where_it_went_boom, shoot_fiting_point
            );

            ode_phiScal_I_LivePlot_append_solver( guess_to_integrate );

            exit(123);
        }

        dvector_copy_to_index_opps(
          guess_to_integrate->y,
          tmp_f1,
          guess_left_n,
          guess_left_indexes
        );

        dvector_copy_to_index_opps(
          guess_to_integrate->y,
          &tmp_f1[guess_left_n],
          guess_right_n,
          guess_right_indexes
        );

        dvector_copy_to_index(
          &v[guess_left_n],
          guess_to_integrate->y,
          guess_right_n,
          guess_right_indexes
        );

        dvector_copy_to_index(
          known_right_values,
          guess_to_integrate->y,
          guess_left_n,
          guess_left_indexes
        );

        guess_to_integrate->x_initial = tmp_x_right;
        guess_to_integrate->x_final = shoot_fiting_point;
        guess_to_integrate->phiScal_inf = shoot_fiting_point;


        // we will try to get rid of any mass or lambda in the scalar field
        // in the right to left integration
        //double \
          //tmp_beta = GV_PARAMETERS_VALUES[1],
          //tmp_m = GV_PARAMETERS_VALUES[2],
          //tmp_lambda = GV_PARAMETERS_VALUES[3];

        //GV_PARAMETERS_VALUES[1] = 0;
        //GV_PARAMETERS_VALUES[2] = 0;
        //GV_PARAMETERS_VALUES[3] = 0;

        ode_phiScal_I_integrate(guess_to_integrate);

        //GV_PARAMETERS_VALUES[1] = tmp_beta;
        //GV_PARAMETERS_VALUES[2] = tmp_m;
        //GV_PARAMETERS_VALUES[3] = tmp_lambda;

        // append only if we are not calculating fadj in the solver
        if(!GV_fadj_switch){
            ode_phiScal_I_LivePlot_append_solver( guess_to_integrate );
        }

        if(guess_to_integrate->did_it_go_boom){

            printf(
              "\n\t\t\t\t\t it whent booommmmm right to left at %e when using fitt %e \n",
              guess_to_integrate->where_it_went_boom, shoot_fiting_point
            );

            ode_phiScal_I_LivePlot_append_solver( guess_to_integrate );

            exit(123);
        }

        dvector_copy_to_index_opps(
          guess_to_integrate->y,
          tmp_f2,
          guess_left_n,
          guess_left_indexes
        );

        dvector_copy_to_index_opps(
          guess_to_integrate->y,
          &tmp_f2[guess_left_n],
          guess_right_n,
          guess_right_indexes
        );

        for(int i=1; i <= n; i++){
            f[i] = tmp_f1[i] - tmp_f2[i];
        }

        iterate = 0;
        guess_to_integrate->phiScal_inf = shoot_fiting_point;
    }

    dvector_copy(tmp_y, guess_to_integrate->y ,guess_to_integrate->eqs_count);

    guess_to_integrate->x_initial = tmp_x_left;
    guess_to_integrate->x_final = tmp_x_right;

    free(tmp_f1);
    free(tmp_f2);
    free(tmp_y);

    return;
}

static void ode_phiScal_I_shoot_fitting_execute(
  ODEsystemStruct *arg, ShootingVarsStruct *shoot_vars, int *newt_check
){

    // setting the variables we will use in the function
    // which calculates the discrepancy caused by initial values
    guess_left_n = shoot_vars->UNknown_left_n;
    guess_left_indexes = shoot_vars->UNknown_left_indexes;
    known_right_values = shoot_vars->known_right_values;

    guess_right_n = shoot_vars->UNknown_right_n;
    guess_right_indexes = shoot_vars->UNknown_right_indexes;
    known_left_values = shoot_vars->known_left_values;

    guess_to_integrate = arg;

    shoot_vars->newt_n = shoot_vars->UNknown_left_n + shoot_vars->UNknown_right_n;

    dvector_copy(
      shoot_vars->UNknown_left_values,
      shoot_vars->newt_v,
      shoot_vars->UNknown_left_n
    );

    dvector_copy(
      shoot_vars->UNknown_right_values,
      &shoot_vars->newt_v[shoot_vars->UNknown_left_n],
      shoot_vars->UNknown_right_n
    );

    ode_phiScal_I_LivePlot_open_solver(arg);

    if(shoot_fiting_point > arg->x_final){
        printf(
          "\n ODE_phiScal_I something fisshy fitting greater than end "
          "%e > %e \n",
          shoot_fiting_point, arg->x_final
        );
        exit(123);
    }

    newt(
      shoot_vars->newt_v,
      shoot_vars->newt_n,
      newt_check,
      &ode_phiScal_I_shooting_phiScal_I_fitting
    );

    //// I tried to make regular shoot to fine tune the guessing values
    // but it turns out that the method on its own is good enough
    // the thing which is problematic is when nan or inf is present
    // and will eventually break the methods used

    dvector_copy_to_index(
      shoot_vars->newt_v,
      arg->y,
      guess_left_n,
      guess_left_indexes
    );

    dvector_copy_to_index(
      known_left_values,
      arg->y,
      guess_right_n,
      guess_right_indexes
    );

    //printf("\n\n\n v");
    //for(int i=1; i <= shoot_vars->newt_n; i++){
        //printf(" %e ", shoot_vars->newt_v[i]);
    //}

    //printf("\n y");
    //for(int i=1; i <= arg->eqs_count; i++){
        //printf(" %e ", arg->y[i]);
    //}

    newt(
      shoot_vars->newt_v,
      guess_left_n,
      newt_check,
      &ode_phiScal_I_shooting_phiScal_I_regular
    );

    //printf("\n y");

    //for(int i=1; i <= arg->eqs_count; i++){
        //printf(" %e ", arg->y[i]);
    //}

    dvector_copy_to_index(
      shoot_vars->newt_v,
      arg->y,
      guess_left_n,
      guess_left_indexes
    );

    dvector_copy_to_index(
      known_left_values,
      arg->y,
      guess_right_n,
      guess_right_indexes
    );

    //printf("\n v");
    //for(int i=1; i <= shoot_vars->newt_n; i++){
        //printf(" %e ", shoot_vars->newt_v[i]);
    //}

    ode_phiScal_I_integrate(guess_to_integrate);

    dvector_copy_to_index_opps(
      arg->y,
      &shoot_vars->newt_v[guess_left_n],
      guess_right_n,
      guess_right_indexes
    );

    //printf("\n v");
    //for(int i=1; i <= shoot_vars->newt_n; i++){
        //printf(" %e ", shoot_vars->newt_v[i]);
    //}

    //sleep(1);

    guess_left_indexes = NULL;
    guess_right_indexes = NULL;
    shoot_free_indexes = NULL;
    known_left_values = NULL;
    known_right_values = NULL;
    guess_to_integrate = NULL;

    return;
}

// get the guessed values
// compare to the previous guessed values to evaluate who has the minimal difference
// save the new guessed over the previous so we can make the same next time
static int minmal_difference(int newt_n, double *newt_v, double *tmp_v){

    int diff = -30;

    for(int k = 1;k <= newt_n; k++){

        int tmp_diff = (int)floor( log10( fabs( newt_v[k] - tmp_v[k] + tiny ) ) );

        diff = tmp_diff > diff ? tmp_diff : diff;

        // will do this in the function
        //tmp_v[k] = newt_v[k];
    }
    return diff;
}

static int boundary_check(double *y, ShootingVarsStruct *shoot_vars){

    for(int i = 1; i <= shoot_vars->UNknown_left_n; i++){
        if( fabs(y[shoot_vars->UNknown_left_indexes[i]] - shoot_vars->known_right_values[i]) > 1e-8 ){
            printf("\n %d %e %d %e \n",shoot_vars->UNknown_left_indexes[i], y[shoot_vars->UNknown_left_indexes[i]], i, shoot_vars->known_right_values[i]);
            return 0;
        }
    }

    return 1;
}


static void ode_phiScal_I_change_central_value(
    ODEsystemStruct *arg, ShootingVarsStruct *shoot_vars
){
    const char \
        function_path[] = "ode_phiScal_I.c ode_phiScal_I_change_central_value: ", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_I_change_central_value){
        printf("%s %s starting \n", identation, function_path);
    }

    double \
      max_val = 1,
      min_val = 0,
      current = min_val,
      *tmp_y = dvector(1, arg->eqs_count);

    if(DEBUGGING_ode_phiScal_I_change_central_value){
        printf("%s %s dvector_copy \n", identation, function_path);
    }
    dvector_copy(arg->y, tmp_y ,arg->eqs_count);

    arg->initial_y_current = arg->initial_y_start;

    printf(
        "\t\t\t \n\t beta = %.3e"
        "\t m = %.3e"
        "\t lambda = %.3e \n",
        arg->free_parmeters_values[1],
        arg->free_parmeters_values[2],
        arg->free_parmeters_values[3]
    );

    int iterate = 1, newt_check;
    double \
      phiScal_c = shoot_vars->UNknown_left_values[1],
      dphiScal_c_dr = shoot_vars->UNknown_right_values[1],
      tmp_x_final = ODE_INDEP_FINAL,
      tmp_fitting_point = GV_PARAMETERS_VALUES[2] ? \
            (FITTING_MASSIVE_COEF/GV_PARAMETERS_VALUES[2]) : (arg->x_final - arg->x_initial)/2,
      // copy the initial guesses into this dummy variable
      // it will be used later to check if there is significant change
      // when increasing the infinities
      *tmp_v = dvector(1,shoot_vars->UNknown_left_n + shoot_vars->UNknown_right_n);

    // loop to go over initial pressurises
    while( iterate ){

        printf("\n p_c = %.3e \n", arg->initial_y_current);

        // we will increase the inf for fine tunning the guessed values
        // so if it was icnearsed, we need to take the new value
        // same holds for the fitting point
        arg->x_final = tmp_x_final;
        shoot_fiting_point = tmp_fitting_point;

        dvector_copy(
          shoot_vars->UNknown_left_values,
          tmp_v,
          shoot_vars->UNknown_left_n
        );

        dvector_copy(
          shoot_vars->UNknown_right_values,
          &tmp_v[shoot_vars->UNknown_left_n],
          shoot_vars->UNknown_right_n
        );

        // loop to increase the infinity until significant difference
        // in the previous and guessed values vanishes
        while( 1 ){

            // initialize the ode system and set the current initial pressure
            dvector_copy(tmp_y,arg->y ,arg->eqs_count);
            arg->y[arg->index_of_y_to_change] = arg->initial_y_current;

            if(DEBUGGING_ode_phiScal_I_change_central_value){

                printf("%s %s dvector_copy \n", identation, function_path);

                printf(
                  "%s %s [%d] %s_init = %.3e  \n",
                  identation, function_path, arg->index_of_y_to_change,
                  arg->name_vars[arg->index_of_y_to_change],
                  arg->y[arg->index_of_y_to_change]
                );
            }

            printf(
              "\n\t inf %e; phiScal_guess = %e \n",
              arg->x_final, shoot_vars->UNknown_left_values[1]
            );

            dvector_copy(shoot_vars->newt_v, tmp_v, shoot_vars->newt_n);

            // TODO make sure this is working also fine
            //ode_phiScal_I_shoot_regular_execute( arg, shoot_vars );

            // I prefer fitting for the moment
            ode_phiScal_I_shoot_fitting_execute( arg, shoot_vars, &newt_check);

            printf(
              "\n\t\t --> \t newt check %d; phiScal_v = %e \n",
              newt_check,
              shoot_vars->newt_v[1]
            );

            // at the end of ode_phiScal_I_shoot_fitting_execute
            // we integrate to check if using the guessed values on left side
            // will produce the wanted values on the right
            // if not most probably there was some nan or inf, so we need to reduce
            // the fitting and infinity
            if( !boundary_check(arg->y, shoot_vars) ){
                printf("\n Boundary problem !!! \n");

                shoot_fiting_point = tmp_fitting_point;
                //exit(123);

            // checking if we have significant difference in the guessed values
            // if the least difference is less than the desired MIN_DIFF_POW
            // our work here is done
            // if not increase the infinity with INCREASE_COEF of the current infinity
            // and if fitting is used, increase it using the same method
            }else if(
              minmal_difference( shoot_vars->newt_n, shoot_vars->newt_v, tmp_v ) <= MIN_DIFF_POW
            ){
                dvector_copy(tmp_v, shoot_vars->newt_v, shoot_vars->newt_n);

                arg->x_final = tmp_x_final;
                shoot_fiting_point = tmp_fitting_point;
                break;
            }else{
                // we are not interested in arbitrary small scalar fields
                // small scalar field is such that has initial value of power less MIN_PHISCAL_POW
                if( floor( log10( fabs( shoot_vars->newt_v[1] ) ) ) <= MIN_PHISCAL_POW ){

                    shoot_vars->newt_v[1] = 0;
                    // lets set the derivative of the scalar field to 0
                    shoot_vars->newt_v[4] = 0;

                    // we shall not keep track of the difference in the guessed values
                    tmp_v[1] = 0;
                    tmp_v[4] = 0;
                }

                // if criteria not met we increase a little bit the infinity
                // save the old infinities and increase them
                tmp_x_final = arg->x_final;
                tmp_fitting_point = shoot_fiting_point;

                arg->x_final += arg->x_final*INCREASE_COEF;
                shoot_fiting_point += shoot_fiting_point*INCREASE_COEF;

                dvector_copy(
                  shoot_vars->newt_v,
                  shoot_vars->UNknown_left_values,
                  shoot_vars->UNknown_left_n
                );

                dvector_copy(
                  &shoot_vars->newt_v[shoot_vars->UNknown_left_n],
                  shoot_vars->UNknown_right_values,
                  shoot_vars->UNknown_right_n
                );
            }
        }

        ode_phiScal_I_ResultFile_append_shootfitting( arg, shoot_vars );

        if(fabs(shoot_vars->newt_v[1]) > fabs(phiScal_c)){
            printf(
              "\n\t changing the initial guess from %e to %e \n",
              phiScal_c,
              shoot_vars->newt_v[1]
            );
            phiScal_c = shoot_vars->newt_v[1];
            shoot_vars->UNknown_left_values[1] = shoot_vars->newt_v[1];
            shoot_vars->UNknown_right_values[1] = shoot_vars->newt_v[4];
            dphiScal_c_dr = shoot_vars->newt_v[4];
        }else if(
          floor( log10( fabs( shoot_vars->newt_v[1] + tiny ) ) ) >= -2
        ){
            phiScal_c = shoot_vars->newt_v[1];
            shoot_vars->UNknown_left_values[1] = shoot_vars->newt_v[1];
            shoot_vars->UNknown_right_values[1] = shoot_vars->newt_v[4];
            dphiScal_c_dr = shoot_vars->newt_v[4];

            printf(
              "\n\t reducing from %e to %e \n",
              shoot_vars->UNknown_left_values[1],
              shoot_vars->newt_v[1]
            );
        }else{
            shoot_vars->UNknown_left_values[1] = phiScal_c;
            shoot_vars->UNknown_right_values[1] = dphiScal_c_dr;
        }

        int current_power = (int)floor(log10(arg->initial_y_current));
        switch(current_power){
            case -5:
                arg->initial_y_current += 2*pow10( current_power - 1 );
                break;
            case -4:
                arg->initial_y_current += 1*pow10( current_power - 1 );
                break;
            case -3:
                arg->initial_y_current += 5*pow10( current_power - 2 );
                break;
            case -2:
                arg->initial_y_current += 5*pow10( current_power - 2 );
                break;
            default:
                printf("ode_phiscal_i line 1372 %d not known terminating", current_power);
                exit(1376);
                break;
        }

        //arg->initial_y_current +=
          //arg->initial_y_step*pow10( floor(log10(arg->initial_y_current))-1 ) ;

        if( arg->initial_y_current > arg->initial_y_end ){
            iterate = 0;
        }else{
            iterate = 1;
        }
    }

    arg->x_final = ODE_INDEP_FINAL;
    free(tmp_y);
    free(tmp_v);

    return;
}

void ode_phiScal_I_compute_parameters( ODEsystemStruct *arg ){

    const char \
        function_path[] = "ode_phiScal_I.c ode_phiScal_I_compute_parameters: ", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_I_compute_parameters){
        printf("%s %s starting \n", identation, function_path);
    }

    double *tmp_y;

    tmp_y = dvector(1, arg->eqs_count);

    if(DEBUGGING_ode_phiScal_I_compute_parameters){
        printf("%s %s dvector_copy %d \n", identation, function_path, arg->eqs_count);
    }
    dvector_copy(arg->y, tmp_y ,arg->eqs_count);

    // setting the free parameters of the ODE system
    for(
      int parm_beta = 1;
      parm_beta <= arg->free_parameters_count_each[1];
      parm_beta++
    ){
        for(
          int parm_m = 1;
          parm_m <= arg->free_parameters_count_each[2];
          parm_m++
        ){
            for(
              int parm_lambda = 1;
              parm_lambda <= arg->free_parameters_count_each[3];
              parm_lambda++
            ){
                arg->free_parmeters_values[1] = \
                  arg->free_parmeters_values_all[1][parm_beta];
                arg->free_parmeters_values[2] = \
                  arg->free_parmeters_values_all[2][parm_m];
                arg->free_parmeters_values[3] = \
                  arg->free_parmeters_values_all[3][parm_lambda];

                if(DEBUGGING_ode_phiScal_I_compute_parameters){
                    printf(
                      "%s %s dvector_copy; \n"
                      "\t %s = %.3e \n"
                      "\t %s = %.3e \n"
                      "\t %s = %.3e \n",
                      identation, function_path,
                      arg->free_parmeters_names[1],
                      arg->free_parmeters_values[1],
                      arg->free_parmeters_names[2],
                      arg->free_parmeters_values[2],
                      arg->free_parmeters_names[3],
                      arg->free_parmeters_values[3]
                    );
                }

                dvector_copy(tmp_y,arg->y ,arg->eqs_count);

                eos = calloc(1,sizeof(EOSmodelInfoStruct));
                eos_init(&eos);
                arg->eoseq_name = eos->model_name;

                ShootingVarsStruct *shoot_vars;
                shoot_vars = calloc(1,sizeof(ShootingVarsStruct));
                shooting_phiScal_I_init(&shoot_vars);
                shooting_phiScal_I_check(shoot_vars, arg);

                ode_phiScal_I_ResultFile_open(arg);
                ode_phiScal_I_LivePlot_open_solver(arg);

                ode_phiScal_I_info_print_stdout(arg);
                ode_phiScal_I_info_print_ResultFile( arg );

                odeint_info_print_stdout();
                odeint_info_print_ResultFile(
                  open_file_to_APPEND_ResultFile(arg)
                );

                newt_info_print_stdout();
                newt_info_print_ResultFile(open_file_to_APPEND_ResultFile(arg));

                shooting_phiScal_I_info_print_stdout(shoot_vars);
                shooting_phiScal_I_info_print_ResultFile(
                  shoot_vars,open_file_to_APPEND_ResultFile(arg)
                );

                // after the parameters were set we can start calculating
                // the different stars
                ode_phiScal_I_change_central_value(arg, shoot_vars);

                eos_free(&eos);
                shooting_phiScal_I_free(&shoot_vars);
            }
        }
    }

    free(tmp_y);
    return;
}

#undef DEBUGGING_ode_phiScal_I_foo
#undef DEBUGGING_ode_phiScal_I_init
#undef DEBUGGING_ode_phiScal_I_free
#undef DEBUGGING_ode_phiScal_I_integrate
#undef ODE_VARS_NAME_LENGHT
#undef MAX_ARR_SIZE
#undef ODE_INDEX
#undef ODE_DESCRIPTION
#undef ODE_EQS_COUNT
#undef ODE_Y_INIT_VAL
#undef ODE_NAME_INDEP
#undef ODE_NAME_DEP
#undef ODE_POINTS_COUNT
#undef ODE_ODE_SCALLING
#undef ODE_RKQS_SCALLING
#undef ODE_FREE_PARM_COUNT_ALL
#undef ODE_FREE_PARM_NAME
#undef ODE_FREE_PARM_COUNT_1
#undef ODE_FREE_PARM_VALS_1
#undef ODE_FREE_PARM_COUNT_2
#undef ODE_FREE_PARM_VALS_2
#undef ODE_FREE_PARM_COUNT_3
#undef ODE_FREE_PARM_VALS_3
#undef ODE_INDEP_INIT
#undef ODE_INDEP_FINAL
#undef Y_INDEX_CHANGE
#undef INITIAL_Y_START
#undef INITIAL_Y_END
#undef INITIAL_Y_STEP
#undef DEBUGGING_ode_phiScal_I_compute_parameters
#undef DEBUGGING_ode_phiScal_I_change_central_value
