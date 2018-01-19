#include "ExternalHeaders.h"

#define DEBUGGING_ode_logistics_foo 0
#define DEBUGGING_ode_logistics_init 0
#define DEBUGGING_ode_logistics_free 0
#define DEBUGGING_ode_logistics_integrate 0

#define ODE_VARS_NAME_LENGHT 16

static double *GV_PARAMETERS_VALUES;

void ode_logistics_foo( double x, double *y, double *dydx ){

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
    (*arg)->index = 0;

    if(DEBUGGING_ode_logistics_init){
        printf( \
          "%s %s ODE system index %d \n", \
          identation, function_path, (*arg)->index \
        );
    };

    // total number of ODE equations
    (*arg)->eqs_count = 1;
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
        10.0 // initial population
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

    // variable names per index
    // the 0 index is the name of the independent variable
    const char *ode_names[] = {
        "t", // time
        "N" // population
    };

    (*arg)->name_vars = calloc( (size_t)(*arg)->eqs_count, sizeof(char*));

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
    strcpy((*arg)->name_system, "ode 0, 1d, logistic equation");

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s System name and description: %s \n", \
          identation,function_path, (*arg)->name_system \
        );
    }

    // amount of points we want to save and print, if any
    (*arg)->points_count = 1e2;
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
    (*arg)->odeint_scaling_method = 0;
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
    (*arg)->rkqs_step_method = 0;
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
    (*arg)->free_parmeters_count_all = 1;
    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s We have %d free parameters in the system \n", \
          identation,function_path, (*arg)->free_parmeters_count_all \
        );
    }

    const char *free_parameters_names[] = {
        "1e10", // empty since it will not be red
        "r" // rate of maximum population growth
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
    double free_parmeters_values_all[][20] = {
        { 1e10 },  // index starts from 1 so this will not be taken into account
        { 1e10, -10, -5, 0, 5, 10 }  // value for r, first value does not taken into account
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
        5 // count all values of parameter 1 we want to check
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
    (*arg)->x_initial = 0;

    // until we reach this one
    (*arg)->x_final = 10;

    if(DEBUGGING_ode_logistics_init){
        printf(\
          "%s %s Independent variable %s initial %.3e and final %.3e value \n", \
          identation,function_path, (*arg)->name_vars[0], \
          (*arg)->x_initial, (*arg)->x_final \
        );
    }

    // odeint integration accuracy parameter
    (*arg)->eps = 1e-12;

    // odeint initial step;
    (*arg)->h1 = 1e-30;

    // odeint minimum step (can be zero)
    (*arg)->hmin = 0;

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

    if( (*arg)->eqs_count ){
        free_dvector( (*arg)->y, 1, (*arg)->eqs_count );
    }else{
        free( (*arg)->y );
    }

    if( (*arg)->points_count ){
        free_dvector( (*arg)->points_x, 1, (*arg)->points_count );

        free_dmatrix( \
            (*arg)->points_y, \
            1, (*arg)->eqs_count, \
            1, (*arg)->points_count \
        );
    }else{
        free((*arg)->points_x);
        free((*arg)->points_y);
    }

    free(*arg);

    return;
}

static void ode_logistics_integrate_info_print( ODEsystemStruct *arg ){

    printf(
        "\n System info: \n \
        index: %d \n \
        description: %s \n",
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

    printf("\n Free parameters info: \n");
    for(int i=1; i<=arg->eqs_count; i++){
        printf("\t value for %s = %.3e \n",
          arg->free_parmeters_names[i], arg->free_parmeters_values[i]
        );
    }

    printf(
        "\n Integrator info: \n \
        name: %s \n \
        scaling method: %s \n \
        rkqs step method: %s \n \
        initial step %.3e \n \
        minimal step %3e \n \
        desired accuracy %.3e \n \
        desired points to record %d \n\n",
        \
        "odeint",
        arg->odeint_scaling_method_description,
        arg->rkqs_step_method_description,
        arg->h1, arg->hmin, arg->eps, arg->points_count
    );
}

void ode_logistics_integrate( ODEsystemStruct **arg ){

    const char function_path[] = "ODE_logistic.c ode_logistics_integrate", \
               identation[] = "\n";

    if(DEBUGGING_ode_logistics_integrate){
        printf("%s %s starting \n", identation, function_path);
    }

    double *y;

    y = dvector(1, (*arg)->eqs_count);
    dvector_copy((*arg)->y, y, (*arg)->eqs_count);

    kmax = (*arg)->points_count;
    xp = (*arg)->points_x;
    yp = (*arg)->points_y;

    ode_logistics_integrate_info_print(*arg);

    odeint(\
        y, (*arg)->eqs_count, (*arg)->x_initial, (*arg)->x_final, \
        (*arg)->eps, (*arg)->h1, (*arg)->hmin, \
        &(*arg)->nok, &(*arg)->nbad, \
        (*arg)->foo \
    );

    if((*arg)->points_count){
        for(int i=1; i <= (*arg)->points_count; i++){
            printf("\n %e %e", xp[i], yp[1][i]);
        }
    }

    free_dvector(y, 1, (*arg)->eqs_count);

    return;
}

#undef DEBUGGING_ode_logistics_foo
#undef DEBUGGING_ode_logistics_init
#undef DEBUGGING_ode_logistics_free
#undef DEBUGGING_ode_logistics_integrate

#undef ODE_VARS_NAME_LENGHT
