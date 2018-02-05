#include "ExternalHeaders.h"

#define DEBUGGING_ode_phiScal_foo 0
#define DEBUGGING_ode_phiScal_init 0
#define DEBUGGING_ode_phiScal_free 0
#define DEBUGGING_ode_phiScal_integrate 1
#define DEBUGGING_ode_phiScal_compute_parameters 0
#define DEBUGGING_ode_phiScal_change_central_value 0

// General variables
#define ODE_VARS_NAME_LENGHT 16
#define MAX_ARR_SIZE 16

// ODE system general info
#define ODE_INDEX 1
#define ODE_DESCRIPTION "STT_phiScal"

// ODE count of equations and initial variable per equation
#define ODE_EQS_COUNT 5
#define ODE_Y_INIT_VAL 0, 0, 5e-4, 0, 0

// ODE names/symbols for independent and for each dependent variable
#define ODE_NAME_INDEP "r"
#define ODE_NAME_DEP "phiScal", "Q", "p", "LambdaMetr", "m"

// odeint mount of points we want to print if any
#define ODE_POINTS_COUNT 1e2

// amount of ode free parameters, name/symbol for each one
#define ODE_FREE_PARM_COUNT_ALL 3
#define ODE_FREE_PARM_NAME "beta", "m", "lambda"

// beta
#define ODE_FREE_PARM_COUNT_1 5
#define ODE_FREE_PARM_VALS_1 0, -4, -4.5, -6, -8

// m
#define ODE_FREE_PARM_COUNT_2 1
// 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1
#define ODE_FREE_PARM_VALS_2 0

// lambda
#define ODE_FREE_PARM_COUNT_3 1
// 1e-2, 1e-1, 1, 10, 100, 500
#define ODE_FREE_PARM_VALS_3 0

// the interval for the independent variable
#define ODE_INDEP_INIT 0
#define ODE_INDEP_FINAL 1e2

// if we are to change vary any of the initial values of the ODE
// index, start and end value, also the step to achieve it
#define Y_INDEX_CHANGE 3
#define INITIAL_Y_START 1e-5
#define INITIAL_Y_END 5e-3
#define INITIAL_Y_STEP 1e-5

static double *GV_PARAMETERS_VALUES;

static void ode_phiScal_foo(double x, double *y, double *dydx){

    double \
        beta_phiScal = GV_PARAMETERS_VALUES[1], \
        m_phiScal = GV_PARAMETERS_VALUES[2], \
        lambda_phiScal = GV_PARAMETERS_VALUES[3], \
        r = x, \
        phiScal = y[1], \
        Q = y[2], \
        p = y[3], \
        LambdaMetr = y[4], \
        m = y[5];

    double \
        Vhat = 2*pow(m_phiScal,2)*pow(phiScal,2) + lambda_phiScal*pow(phiScal,4), \
        Vhat_dphiScal = 4*pow(m_phiScal,2)*phiScal + 4*lambda_phiScal*pow(phiScal,3), \
        alpha = beta_phiScal*phiScal, \
        A = exp(1.0/2 * beta_phiScal * pow(phiScal,2)), \
        exp_2LambdaMetr = exp(2*LambdaMetr );

    double \
        step4_A = pow(A,4);

    //TODO no EOSeq yet
    //EOSeq(EOSdata,p);
    double rho = 0;

    double \
      PhiMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*p*exp_2LambdaMetr
          + pow(Q,2)
          - 1.0/2*exp_2LambdaMetr*Vhat
          + pow(1/r,2)*(exp_2LambdaMetr-1)
        ), \
      phiScal_dr = Q, \
      LambdaMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*rho*exp_2LambdaMetr
          + pow(Q,2)
          + 1.0/2*exp_2LambdaMetr*Vhat
          - pow(1/r,2)*(exp_2LambdaMetr-1)
        ), \
      Q_dr = \
        4*GV_PI*alpha*step4_A*(rho - 3*p)*exp_2LambdaMetr \
        + 1.0/4*Vhat_dphiScal*exp_2LambdaMetr \
        - (PhiMetr_dr - LambdaMetr_dr + 2/r )*Q, \
      p_dr = -( rho + p ) * ( PhiMetr_dr + alpha*Q ), \
      m_dr = \
        pow(r,2) * (
          4 * GV_PI * step4_A* rho
          + 1.0/2 * 1/exp_2LambdaMetr * pow(Q,2)
          + 1.0/4 * Vhat
        );

    dydx[1] = phiScal_dr;
    dydx[2] = Q_dr;
    dydx[3] = p_dr;
    dydx[4] = LambdaMetr_dr;
    dydx[5] = m_dr;

    return;
}

void ode_phiScal_init( ODEsystemStruct **arg ){

    const char \
        function_path[] = "ODE_phiScal.c ode_phiScal_init: ", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_init){
        printf("%s %s starting \n", identation, function_path);
    }

    if( ( *arg = calloc(1,sizeof(ODEsystemStruct)) ) == NULL ){
        printf("%s %s error in calloc for ode, line 92 \n", identation, function_path);
        exit(92);
    };

    // ODE system index
    (*arg)->index = ODE_INDEX;

    if(DEBUGGING_ode_phiScal_init){
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
    }else if(DEBUGGING_ode_phiScal_init){
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

    if(DEBUGGING_ode_phiScal_init){
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

        if(DEBUGGING_ode_phiScal_init){

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

        if(DEBUGGING_ode_phiScal_init){
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

        if(DEBUGGING_ode_phiScal_init){
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

    if(DEBUGGING_ode_phiScal_init){
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

        xp = (*arg)->points_x;
        yp = (*arg)->points_y;
    }else{

        xp = NULL;
        yp = NULL;

    }

    if(DEBUGGING_ode_phiScal_init){
        printf(\
          "%s %s We want to print %d points \n", \
          identation,function_path, (*arg)->points_count \
        );
    }

    // count OK steps done by the integrator
    (*arg)->nok = 0;

    if(DEBUGGING_ode_phiScal_init){
        printf(\
          "%s %s We have %d OK points \n", \
          identation,function_path, (*arg)->nok \
        );
    }

    // count BAD steps done by the integrator
    (*arg)->nbad = 0;

    if(DEBUGGING_ode_phiScal_init){
        printf(\
          "%s %s We have %d BAD points \n", \
          identation,function_path, (*arg)->nok \
        );
    }

    // how many free parameters the system has
    (*arg)->free_parmeters_count_all = ODE_FREE_PARM_COUNT_ALL;
    if(DEBUGGING_ode_phiScal_init){
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

        if(DEBUGGING_ode_phiScal_init){
            printf(\
              "%s %s [%d] free parameter name is %s\n", \
              identation,function_path, i, (*arg)->free_parmeters_names[i]\
            );
        }
    }

    (*arg)->free_parmeters_values = dvector(1, (*arg)->free_parmeters_count_all);
    GV_PARAMETERS_VALUES = (*arg)->free_parmeters_values;

    if(DEBUGGING_ode_phiScal_init){
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

        if(DEBUGGING_ode_phiScal_init){
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

            if(DEBUGGING_ode_phiScal_init){
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

    if(DEBUGGING_ode_phiScal_init){
        printf(\
          "%s %s Independent variable %s initial %.3e and final %.3e value \n", \
          identation,function_path, (*arg)->name_vars[0], \
          (*arg)->x_initial, (*arg)->x_final \
        );
    }

    (*arg)->foo = &ode_phiScal_foo;

    if(DEBUGGING_ode_phiScal_init){
        printf("%s %s ending \n", identation, function_path);
    }

    return;
}

static void ode_phiScal_info_print_stdout( ODEsystemStruct *arg ){

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
        "\t\t final value: %.3e \n",
        arg->name_vars[0], arg->x_initial, arg->x_final
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

static void ode_phiScal_info_print_ResultFile( ODEsystemStruct *arg, FILE *fp ){

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
        "\t\t final value: %.3e \n",
        arg->name_vars[0], arg->x_initial, arg->x_final
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

static void ode_phiScal_LivePlot_open(ODEsystemStruct *arg){

    FILE *fp = open_file_to_WRITE_LivePlot(arg);

    fclose(fp);

    return;
}

static void ode_phiScal_ResultFile_open(ODEsystemStruct *arg){

    FILE *fp = open_file_to_WRITE_ResultFile(arg);

    fclose(fp);

    return;
}

static void ode_phiScal_LivePlot_append(ODEsystemStruct *arg){

    FILE *fp = open_file_to_APPEND_LivePlot(arg);

    fprintf(fp,"# ");

    for(int i=1; i <= arg->free_parmeters_count_all; i++){
        fprintf(
            fp,
            "y_c = %.3e",
            arg->initial_y_current
        );
    }

    if(arg->points_count){
        for(int i=1; i <= arg->points_count && xp[i]; i++){
            fprintf(fp,"\n%e, %e", xp[i], yp[1][i]);
        }
    }

    fprintf(fp,"\n");

    fclose(fp);

    return;
}

static void ode_phiScal_ResultFile_append(ODEsystemStruct *arg){

    FILE *fp = open_file_to_APPEND_ResultFile(arg);

    fprintf(
      fp,
      "\n%e, %e",
      arg->initial_y_current, arg->y[arg->index_of_y_to_change]
    );

    fclose(fp);

    return;
}

void ode_phiScal_free( ODEsystemStruct **arg ){

    const char \
        function_path[] = "ODE_logistic.c ode_logistics_free", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_free){
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

    free(*arg);

    return;
}

static void ode_phiScal_integrate( ODEsystemStruct *arg ){

    const char function_path[] = "ODE_phiScal.c ode_phiScal_integrate", \
               identation[] = "\n";

    if(DEBUGGING_ode_phiScal_integrate){
        printf("%s %s starting \n", identation, function_path);
    }

    if(kmax){
        for(int i=1; i <= arg->points_count; i++){
            xp[i] = 0;
            for(int l=1; l <= arg->eqs_count; l++){
                yp[l][i] = 0;
            }
        }
    }else{
        xp = NULL;
        yp = NULL;
    }

    arg->nok = arg->nbad = 0;

    if(DEBUGGING_ode_phiScal_integrate){
        printf("\n");
        for(int i=1; i<=arg->eqs_count;i++){
            printf(" %s_init = %.3e ", arg->name_vars[i], arg->y[i]);
        }
    }

    odeint(
      arg->y, arg->eqs_count,
      arg->x_initial, arg->x_final,
      &arg->nok, &arg->nbad,
      arg->foo
    );

    if(DEBUGGING_ode_phiScal_integrate){
        printf("\n");
        for(int i=1; i<=arg->eqs_count;i++){
            printf(" %s_init = %.3e ", arg->name_vars[i], arg->y[i]);
        }
    }

    return;
}

static void ode_phiScal_change_central_value(
    ODEsystemStruct *arg, ShootingVarsStruct *shoot_regular_vars
){
    const char \
        function_path[] = "ODE_phiScal.c ode_phiScal_change_central_value: ", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_change_central_value){
        printf("%s %s starting \n", identation, function_path);
    }

    double max_val = 1, min_val = 0, current = min_val, step = 0.05, *tmp_y;

    tmp_y = dvector(1, arg->eqs_count);

    if(DEBUGGING_ode_phiScal_change_central_value){
        printf("%s %s dvector_copy \n", identation, function_path);
    }
    dvector_copy(arg->y, tmp_y ,arg->eqs_count);

    arg->initial_y_current = arg->initial_y_start;

    printf(
        "\n beta = %.3e"
        "\t m = %.3e"
        "\t lambda = %.3e",
        arg->free_parmeters_values[1],
        arg->free_parmeters_values[2],
        arg->free_parmeters_values[3]
    );

    int iterate = 1;
    while( iterate ){

        if(DEBUGGING_ode_phiScal_change_central_value){
            printf("%s %s dvector_copy \n", identation, function_path);
        }

        dvector_copy(tmp_y,arg->y ,arg->eqs_count);

        arg->y[arg->index_of_y_to_change] = arg->initial_y_current;

        if(DEBUGGING_ode_phiScal_change_central_value){
            printf(
              "%s %s [%d] %s_init = %.3e  \n",
              identation, function_path, arg->index_of_y_to_change,
              arg->name_vars[arg->index_of_y_to_change],
              arg->y[arg->index_of_y_to_change]
            );
        }

        ode_phiScal_integrate(arg);

        printf(
          "\n p_c = %.3e; \t phi_c = %.3e -> phi_inf = %.3e",
          arg->initial_y_current,
          shoot_regular_vars->newt_v[1],
          arg->y[1]
        );

        if( arg->initial_y_current <= arg->initial_y_end ){
            arg->initial_y_current += arg->initial_y_step;
            iterate = 1;
        }else{
            iterate = 0;
        }

        unsigned short time_interval = 1;
        printf(
          "\n\t\t\t small time sleep %d s\n",
          time_interval
        );
        sleep(time_interval);
    }

    free(tmp_y);

    return;
}

void ode_phiScal_compute_parameters( ODEsystemStruct *arg ){

    const char \
        function_path[] = "ODE_phiScal.c ode_phiScal_compute_parameters: ", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_phiScal_compute_parameters){
        printf("%s %s starting \n", identation, function_path);
    }

    double *tmp_y;

    tmp_y = dvector(1, arg->eqs_count);

    if(DEBUGGING_ode_phiScal_compute_parameters){
        printf("%s %s dvector_copy %d \n", identation, function_path, arg->eqs_count);
    }
    dvector_copy(arg->y, tmp_y ,arg->eqs_count);

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

                if(DEBUGGING_ode_phiScal_compute_parameters){
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

                ShootingVarsStruct *shoot_regular_vars;
                shoot_regular_vars = calloc(1,sizeof(ShootingVarsStruct));
                shooting_regular_init(&shoot_regular_vars);
                shooting_regular_check(shoot_regular_vars, arg);

                ode_phiScal_LivePlot_open(arg);
                ode_phiScal_ResultFile_open(arg);

                ode_phiScal_info_print_stdout(arg);
                ode_phiScal_info_print_ResultFile(
                  arg,
                  open_file_to_APPEND_ResultFile(arg)
                );

                odeint_info_print_stdout();
                odeint_info_print_ResultFile(
                  open_file_to_APPEND_ResultFile(arg)
                );

                newt_info_print_stdout();
                newt_info_print_ResultFile(open_file_to_APPEND_ResultFile(arg));

                shooting_regular_info_print_stdout(shoot_regular_vars);
                shooting_regular_info_print_ResultFile(
                  shoot_regular_vars,open_file_to_APPEND_ResultFile(arg)
                );

                ode_phiScal_change_central_value(arg, shoot_regular_vars);

                unsigned short time_interval = 10;
                printf(
                  "\n\n TIME TO SLEEP FOR RESULT PLOT CHECK for %d s \n\n",
                  time_interval
                );
                sleep(time_interval);
            }
        }
    }

    free(tmp_y);
    return;
}

#undef DEBUGGING_ode_phiScal_foo
#undef DEBUGGING_ode_phiScal_init
#undef DEBUGGING_ode_phiScal_free
#undef DEBUGGING_ode_phiScal_integrate
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
#undef DEBUGGING_ode_phiScal_compute_parameters
#undef DEBUGGING_ode_phiScal_change_central_value
