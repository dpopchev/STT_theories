#include "ExternalHeaders.h"

#define DEBUGGING_ode_logistics_foo 0
#define DEBUGGING_ode_logistics_init 0
#define DEBUGGING_ode_logistics_free 0
#define DEBUGGING_ode_logistics_integrate 0

// General variables
#define ODE_VARS_NAME_LENGHT 16
#define MAX_ARR_SIZE 16

// ODE system general info
#define ODE_INDEX 0
#define ODE_DESCRIPTION "LogisticEq"

// ODE count of equations and initial variable per equation
#define ODE_EQS_COUNT 1
#define ODE_Y_INIT_VAL_1 1

// ODE names/symbols for independent and for each dependent variable
#define ODE_NAME_INDEP "t"
#define ODE_NAME_DEP_1 "N"

// odeint mount of points we want to print if any
#define ODE_POINTS_COUNT 1e2

// odeint will use those scaling and rkqs step size methods
#define ODE_ODE_SCALLING 0
#define ODE_RKQS_SCALLING 0

// amount of ode free parameters, name/symbol for each one
#define ODE_FREE_PARM_COUNT_ALL 1
#define ODE_FREE_PARM_NAME_1 "r"

// if we will search for more values for some parameters, then how many
// and what will be the values be
#define ODE_FREE_PARM_COUNT_1 5
#define ODE_FREE_PARM_VALS_1 -1, -0.5, 0, 0.5, 1

// the interval for the independent variable
#define ODE_INDEP_INIT 0
#define ODE_INDEP_FINAL 10

// odeint mandatory parameters for accuracy, initial step and min step size
#define ODE_INTEGR_EPS 1e-12
#define ODE_INTEGR_H1 1e-30
#define ODE_INTEGR_HMIN 0

#define Y_INDEX_CHANGE 1
#define INITIAL_Y_START 0
#define INITIAL_Y_END 1
#define INITIAL_Y_STEP 0.05

static double *GV_PARAMETERS_VALUES;

static void ode_logistics_foo( double x, double *y, double *dydx ){

    const char \
        function_path[] = "All_ODE.c ode_logistics_init", \
        identation[] = "\n\t";

    double \
        // free parameter for the logistic equation
        r = GV_PARAMETERS_VALUES[1], \
        // time variable
        t = x, \
        // population
        N = y[1], \
        // result variable for this iteration
        result = 0.0;

    dydx[1] = r*N*(1-N);

    return;
}

void ode_logistics_init( ODEsystemStruct **arg ){

    const char \
        function_path[] = "All_ODE.c ode_logistics_init: ", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_logistics_init){
        printf("%s %s starting \n", identation, function_path);
    }

    if( ( *arg = calloc(1,sizeof(ODEsystemStruct)) ) == NULL ){
        printf("%s %s error in calloc for ode \n", identation, function_path);
        exit(123);
    };

    // ODE system index
    (*arg)->index = ODE_INDEX;

    if(DEBUGGING_ode_logistics_init){
        printf( \
          "%s %s ODE system index %d \n", \
          identation, function_path, (*arg)->index \
        );
    };

    // total number of ODE equations
    (*arg)->eqs_count = ODE_EQS_COUNT;
    if(!(*arg)->eqs_count){
        printf(\
          "%s %s HOW COULD IT BE ODE SYSTEM WITH SIZE 0!!!!", \
          identation,function_path \
        );

        exit(123);
    }else if(DEBUGGING_ode_logistics_init){
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
        ODE_Y_INIT_VAL_1 // initial population
    };

    (*arg)->y = dvector(1, (*arg)->eqs_count);
    dvector_copy(y, (*arg)->y, (*arg)->eqs_count);

    if(DEBUGGING_ode_logistics_init){
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

        if(DEBUGGING_ode_logistics_init){

            printf(
              "%s %s We will change initial value of y[%d] from %.3e to %.3e \n",
              identation,function_path, (*arg)->index_of_y_to_change,
              (*arg)->initial_y_start, (*arg)->initial_y_end
            );
        }
    }else{
        (*arg)->initial_y_start = 0;
        (*arg)->initial_y_end = 0;
        (*arg)->initial_y_current = 0;

        if(DEBUGGING_ode_logistics_init){
            printf("%s %s We will NOT change \n",identation,function_path);
        }
    }

    // variable names per index
    // the 0 index is the name of the independent variable
    const char *ode_names[] = {
        ODE_NAME_INDEP, // time
        ODE_NAME_DEP_1 // population
    };

    (*arg)->name_vars = calloc( (size_t)(*arg)->eqs_count+1, sizeof(char*));

    for(int i=0; i <=(*arg)->eqs_count; i++ ){
        (*arg)->name_vars[i] = calloc( ODE_VARS_NAME_LENGHT, sizeof(char));

        strcpy((*arg)->name_vars[i],ode_names[i]);

        if(DEBUGGING_ode_logistics_init){
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

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s System name and description: %s \n", \
          identation,function_path, (*arg)->name_system \
        );
    }

    // amount of points we want to save and print, if any
    (*arg)->points_count = ODE_POINTS_COUNT;
    if( (*arg)->points_count){
        (*arg)->points_x = dvector(1, (*arg)->points_count);
        (*arg)->points_y = dmatrix(1, (*arg)->eqs_count, 1, (*arg)->points_count);
    }

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s We want to print %d points \n", \
          identation,function_path, (*arg)->points_count \
        );
    }

    // count OK steps done by the integrator
    (*arg)->nok = 0;

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s We have %d OK points \n", \
          identation,function_path, (*arg)->nok \
        );
    }

    // count BAD steps done by the integrator
    (*arg)->nbad = 0;

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s We have %d BAD points \n", \
          identation,function_path, (*arg)->nok \
        );
    }

    // scaling method (odeint), see Integrator_odeint.c for more details
    // default is 0
    (*arg)->odeint_scaling_method = ODE_ODE_SCALLING;
    strcpy(\
        (*arg)->odeint_scaling_method_description, \
        ODEINT_SCALING_METHOD_DESCRIPTION[(*arg)->odeint_scaling_method] \
    );

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s We choose %d scalling method for ode int: %s \n", \
          identation,function_path, \
          (*arg)->odeint_scaling_method, (*arg)->odeint_scaling_method_description \
        );
    }

    // method to evaluate the new value for x in rkqs
    // see odeint.c for more details
    // default is 0
    (*arg)->rkqs_step_method = ODE_RKQS_SCALLING;
    strcpy(\
        (*arg)->rkqs_step_method_description, \
        RKQS_STEP_METHOD_DESCRIPTION[(*arg)->rkqs_step_method] \
    );

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s We choose %d scalling method for rqs: %s \n", \
          identation,function_path, \
          (*arg)->rkqs_step_method, (*arg)->rkqs_step_method_description \
        );
    }

    // how many free parameters the system has
    (*arg)->free_parmeters_count_all = ODE_FREE_PARM_COUNT_ALL;
    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s We have %d free parameters in the system \n", \
          identation,function_path, (*arg)->free_parmeters_count_all \
        );
    }

    const char *free_parameters_names[] = {
        "1e10", // empty since it will not be red
        ODE_FREE_PARM_NAME_1 // rate of maximum population growth
    };

    (*arg)->free_parmeters_names = \
        calloc( (size_t)(*arg)->free_parmeters_count_all + 1, sizeof(char*) );

    for(int i=1; i <=(*arg)->free_parmeters_count_all; i++ ){

        (*arg)->free_parmeters_names[i] = calloc( ODE_VARS_NAME_LENGHT, sizeof(char));

        strcpy((*arg)->free_parmeters_names[i],free_parameters_names[i]);

        if(DEBUGGING_ode_logistics_init){
            printf(\
              "%s %s [%d] free parameter name is %s\n", \
              identation,function_path, i, (*arg)->free_parmeters_names[i]\
            );
        }
    }

    (*arg)->free_parmeters_values = dvector(1, (*arg)->free_parmeters_count_all);
    GV_PARAMETERS_VALUES = (*arg)->free_parmeters_values;

    if(DEBUGGING_ode_logistics_init){
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
        { 1e10, ODE_FREE_PARM_VALS_1 }  // value for r, first value does not taken into account
    };

    for(int i=1; i <= (*arg)->free_parmeters_count_all; i++){
        (*arg)->free_parmeters_values[i] = free_parmeters_values_all[i][1];

        if(DEBUGGING_ode_logistics_init){
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
        ODE_FREE_PARM_COUNT_1 // count all values of parameter 1 we want to check
    };

    for(int i=1; i <= (*arg)->free_parmeters_count_all; i++){

        (*arg)->free_parameters_count_each[i] = \
            free_parameters_count_each[i];

        (*arg)->free_parmeters_values_all[i] = \
            calloc( (size_t)(*arg)->free_parameters_count_each[i]+1, sizeof(double));

        for(int l=1; l <= (*arg)->free_parameters_count_each[i]; l++){

            (*arg)->free_parmeters_values_all[i][l] = \
                free_parmeters_values_all[i][l];

            if(DEBUGGING_ode_logistics_init){
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

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s Independent variable %s initial %.3e and final %.3e value \n", \
          identation,function_path, (*arg)->name_vars[0], \
          (*arg)->x_initial, (*arg)->x_final \
        );
    }

    // odeint integration accuracy parameter
    (*arg)->eps = ODE_INTEGR_EPS;

    // odeint initial step;
    (*arg)->h1 = ODE_INTEGR_H1;

    // odeint minimum step (can be zero)
    (*arg)->hmin = ODE_INTEGR_HMIN;

    if(DEBUGGING_ode_logistics_init){
        printf(
          "%s %s Ode integrator epsilon = %.3e, initial step = %.3e, min step = %.3e \n",
          identation,function_path,
          (*arg)->eps, (*arg)->h1, (*arg)->hmin
        );
    }

    (*arg)->foo = &ode_logistics_foo;

    if(DEBUGGING_ode_logistics_init){
        printf("%s %s ending \n", identation, function_path);
    }

    return;
}

void ode_logistics_free( ODEsystemStruct **arg ){

    const char \
        function_path[] = "All_ODE.c ode_logistics_free", \
        identation[] = "\n\t";

    if(DEBUGGING_ode_logistics_free){
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

static void ode_logistics_integrate_info_print_stdout( ODEsystemStruct *arg ){

    printf(
        "\n\n System info: \n \
        index: %d \n \
        name: %s \n",
        arg->index,
        arg->name_system
    );

    printf(
        "\n Independent variable info: \n \
        name: %s \n \
        initial value: %.3e \n \
        final value: %.3e \n",
        arg->name_vars[0], arg->x_initial, arg->x_final
    );

    printf(
        "\n Dependent variables info: \n"
    );
    for(int i=1; i<=arg->eqs_count; i++){
        printf(
            "\t initial value for %s = %.3e \n",
            arg->name_vars[i], arg->y[i]
        );
    }

    if(arg->index_of_y_to_change){
        printf(
          "\t Will change %s from %.3e to %.3e with step %.3e \n",
          arg->name_vars[arg->index_of_y_to_change],
          arg->initial_y_start,
          arg->initial_y_end,
          arg->initial_y_step
        );
    }

    printf("\n Free parameters info: \n");
    for(int i=1; i<=arg->eqs_count; i++){

        printf(
            "\t current value for %s = %.3e \n",
            arg->free_parmeters_names[i], arg->free_parmeters_values[i]
        );

        printf("\t rest values for %s = ", arg->free_parmeters_names[i] );

        for(int l=1; l<=arg->free_parameters_count_each[i]; l++){
            printf(" %.3e ", arg->free_parmeters_values_all[i][l] );
        }
    }

    printf(
        "\n\n Integrator info: \n \
        name: %s \n \
        scaling method: %s \n \
        rkqs step method: %s \n \
        initial step %.3e \n \
        minimal step %.3e \n \
        desired accuracy %.3e \n \
        desired points to record %d \n\n",
        \
        "odeint",
        arg->odeint_scaling_method_description,
        arg->rkqs_step_method_description,
        arg->h1, arg->hmin, arg->eps, arg->points_count
    );

    return;
}

static void ode_logistics_integrate_info_print_ResultFile( ODEsystemStruct *arg ){

    FILE *fp = open_file_to_WRITE_ResultFile(arg);

    fprintf(
        fp,
        "\n System info: \n \
        index: %d \n \
        name: %s \n",
        arg->index,
        arg->name_system
    );

    fprintf(
        fp,
        "\n Independent variable info: \n \
        name: %s \n \
        initial value: %.3e \n \
        final value: %.3e \n",
        arg->name_vars[0], arg->x_initial, arg->x_final
    );

    fprintf(
        fp,
        "\n Dependent variables info: \n"
    );
    for(int i=1; i<=arg->eqs_count; i++){
        fprintf(
            fp,
            "\t initial value for %s = %.3e \n",
            arg->name_vars[i], arg->y[i]
        );
    }

    if(arg->index_of_y_to_change){
        fprintf(
          fp,
          "\t Will change %s from %.3e to %.3e with step %.3e \n",
          arg->name_vars[arg->index_of_y_to_change],
          arg->initial_y_start,
          arg->initial_y_end,
          arg->initial_y_step
        );
    }

    fprintf( fp,"\n Free parameters info: \n");
    for(int i=1; i<=arg->eqs_count; i++){

        fprintf(
            fp,
            "\t current value for %s = %.3e \n",
            arg->free_parmeters_names[i], arg->free_parmeters_values[i]
        );

        fprintf(
            fp,
            "\t rest values for %s = ", arg->free_parmeters_names[i]
        );

        for(int l=1; l<=arg->free_parameters_count_each[i]; l++){
            fprintf( fp," %.3e ", arg->free_parmeters_values_all[i][l] );
        }
    }

    fprintf(
        fp,
        "\n\n Integrator info: \n \
        name: %s \n \
        scaling method: %s \n \
        rkqs step method: %s \n \
        initial step %.3e \n \
        minimal step %.3e \n \
        desired accuracy %.3e \n \
        desired points to record %d \n\n",
        \
        "odeint",
        arg->odeint_scaling_method_description,
        arg->rkqs_step_method_description,
        arg->h1, arg->hmin, arg->eps, arg->points_count
    );

    fprintf(fp, "\n\nN_init, N_final");

    fclose(fp);

    return;
}

static void ode_logistics_integrate( ODEsystemStruct **arg ){

    const char function_path[] = "ODE_logistic.c ode_logistics_integrate", \
               identation[] = "\n";

    if(DEBUGGING_ode_logistics_integrate){
        printf("%s %s starting \n", identation, function_path);
    }

    double *y;

    y = dvector(1, (*arg)->eqs_count);
    dvector_copy((*arg)->y, y, (*arg)->eqs_count);

    kmax = (*arg)->points_count;

    if(kmax){

        xp = (*arg)->points_x;
        yp = (*arg)->points_y;

        for(int i=1; i <= (*arg)->points_count; i++){
            xp[i] = 0;
            for(int l=1; l <= (*arg)->eqs_count; l++){
                yp[l][i] = 0;
            }
        }
    }else{
        xp = NULL;
        yp = NULL;
    }

    //ode_logistics_integrate_info_print(*arg);

    odeint(
      y, (*arg)->eqs_count, (*arg)->x_initial, (*arg)->x_final,
      (*arg)->eps, (*arg)->h1, (*arg)->hmin,
      &(*arg)->nok, &(*arg)->nbad,
      (*arg)->foo
    );

    //if((*arg)->points_count){
        //for(int i=1; i <= (*arg)->points_count; i++){
            //printf("\n %e %e", xp[i], yp[1][i]);
        //}
    //}

    dvector_copy(y, (*arg)->y, (*arg)->eqs_count);

    free_dvector(y, 1, (*arg)->eqs_count);

    return;
}

static void ode_logistics_LivePlot_open(ODEsystemStruct *arg ){

    FILE *fp = open_file_to_WRITE_LivePlot(arg);

    fclose(fp);

    return;
}

static void ode_logistics_LivePlot_append(ODEsystemStruct *arg){

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

static void ode_logistics_ResultFile_append(ODEsystemStruct *arg){

    FILE *fp = open_file_to_APPEND_ResultFile(arg);

    fprintf(
      fp,
      "\n%e, %e",
      arg->initial_y_current, arg->y[arg->index_of_y_to_change]
    );

    fclose(fp);

    return;
}

static void ode_logistics_change_central_value(ODEsystemStruct **arg){

    double max_val = 1, min_val = 0, current = min_val, step = 0.05, *tmp_y;

    tmp_y = dvector(1, (*arg)->eqs_count);
    dvector_copy((*arg)->y, tmp_y ,(*arg)->eqs_count);

    ode_logistics_integrate_info_print_stdout(*arg);

    ode_logistics_LivePlot_open(*arg);

    (*arg)->initial_y_current = (*arg)->initial_y_start;

    int iterate = 1;
    while( iterate ){

        dvector_copy(tmp_y,(*arg)->y ,(*arg)->eqs_count);

        (*arg)->y[(*arg)->index_of_y_to_change] = (*arg)->initial_y_current;

        ode_logistics_integrate(arg);

        printf(
          "\n N_init = %.3e, N_fin = %.3e",
          (*arg)->initial_y_current, (*arg)->y[(*arg)->index_of_y_to_change]
        );

        ode_logistics_LivePlot_append(*arg);

        ode_logistics_ResultFile_append(*arg);

        if( (*arg)->initial_y_current <= (*arg)->initial_y_end ){
            (*arg)->initial_y_current += (*arg)->initial_y_step;
            iterate = 1;
        }else{
            iterate = 0;
        }
    }

    free(tmp_y);
}

void ode_logistics_compute_parameters( ODEsystemStruct **arg ){

    double *tmp_y;

    tmp_y = dvector(1, (*arg)->eqs_count);
    dvector_copy((*arg)->y, tmp_y ,(*arg)->eqs_count);

    for(
      int parm_r = 1;
      parm_r <= (*arg)->free_parameters_count_each[1];
      parm_r++
    ){

        (*arg)->free_parmeters_values[1] = \
          (*arg)->free_parmeters_values_all[1][parm_r];

        dvector_copy(tmp_y,(*arg)->y ,(*arg)->eqs_count);

        //ode_logistics_integrate_info_print_stdout(*arg);
        ode_logistics_integrate_info_print_ResultFile(*arg);

        ode_logistics_LivePlot_open(*arg);

        //ode_logistics_integrate(arg);
        ode_logistics_change_central_value(arg);

        int time_interval = 30;
        printf("\n\n TIME TO SLEEP FOR RESULT CHECK for %d s \n\n",time_interval );
        sleep(time_interval);
    }

    free(tmp_y);
    return;
}

#undef DEBUGGING_ode_logistics_foo
#undef DEBUGGING_ode_logistics_init
#undef DEBUGGING_ode_logistics_free
#undef DEBUGGING_ode_logistics_integrate
#undef ODE_VARS_NAME_LENGHT
#undef MAX_ARR_SIZE
#undef ODE_EQS_COUNT
#undef ODE_Y_INIT_VAL_1
#undef ODE_NAME_INDEP
#undef ODE_NAME_DEP_1
#undef ODE_POINTS_COUNT
#undef ODE_ODE_SCALLING
#undef ODE_RKQS_SCALLING
#undef ODE_FREE_PARM_COUNT_ALL
#undef ODE_FREE_PARM_NAME_1
#undef ODE_FREE_PARM_COUNT_1
#undef ODE_FREE_PARM_VALS_1
#undef ODE_INDEP_INIT
#undef ODE_INDEP_FINAL
#undef ODE_INTEGR_EPS
#undef ODE_INTEGR_H1
#undef ODE_INTEGR_HMIN
#undef Y_INDEX_CHANGE
#undef INITIAL_Y_START
#undef INITIAL_Y_END
#undef INITIAL_Y_STEP
