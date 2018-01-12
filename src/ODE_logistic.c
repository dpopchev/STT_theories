#include "ExternalHeaders.h"

#define DEBUGGING 0

void ode_logistics_foo( double x, double *y, double *dydx ){

    const char \
        function_path[] = "All_ODE.c ode_logistics_init", \
        identation[] = "\n\t";

    double \
        // free parameter for the logistic equation
        r = 5, \
        // time variable
        t = x, \
        // population
        N = y[1], \
        // result variable for this iteration
        result = 0.0;

    result = r*N*(1-N);

    dydx[1] = result;

    return;
}

void ode_logistics_init( ODEsystemStruct **arg ){

    const char \
        function_path[] = "All_ODE.c ode_logistics_init", \
        identation[] = "\n\t";

    if(DEBUGGING){
        printf("%s %s starting \n", identation, function_path);
    }

    if( ( *arg = calloc(1,sizeof(ODEsystemStruct)) ) == NULL ){
        printf("%s %s error in calloc for ode \n", identation, function_path);
        exit(123);
    };

    // ODE system number
    (*arg)->index = 0;

    // total number of ODE equations
    (*arg)->eqs_count = 1;

    // amount of points we want to save and print, if any
    (*arg)->points_count = 1e2;

    if( (*arg)->points_count){
        (*arg)->points_x = dvector(1, (*arg)->points_count);
        (*arg)->points_y = dmatrix(1, (*arg)->eqs_count, 1, (*arg)->points_count);
    }

    // count OK steps done by the integrator
    (*arg)->nok = 0;

    // count BAD steps done by the integrator
    (*arg)->nbad = 0;

    // scaling method (odeint), see Integrator_odeint.c for more details
    // default is 0
    (*arg)->odeint_scaling_method = 0;
    strcpy(\
        (*arg)->odeint_scaling_method_description, \
        ODEINT_SCALING_METHOD_DESCRIPTION[(*arg)->odeint_scaling_method] \
    );

    // method to evaluate the new value for x in rkqs
    // see odeint.c for more details
    // default is 0
    (*arg)->rkqs_step_method = 0;
    strcpy(\
        (*arg)->rkqs_step_method_description, \
        RKQS_STEP_METHOD_DESCRIPTION[(*arg)->rkqs_step_method] \
    );

    // we will integrate starting from this point
    (*arg)->x_initial = 0;

    // until we reach this one
    (*arg)->x_final = 10;

    // odeint integration accuracy parameter
    (*arg)->eps = 1e-12;

    // odeint initial step;
    (*arg)->h1 = 1e-30;

    // odeint minimum step (can be zero)
    (*arg)->hmin = 0;

    // initial values for *y
    // they should be as many as (*arg)->eqs_count
    double y[] = { \
        1e10, // THIS VALUE DOES NOT TAKE ANY EFFECT, since dvector starts from index 1
        10.0 \
    };

    (*arg)->y = dvector(1, (*arg)->eqs_count);
    dvector_copy(y, (*arg)->y, (*arg)->eqs_count);

    strcpy((*arg)->name_system, "ode 0, 1d, logistic equation");
    strcpy((*arg)->name_vars, "N");

    if(DEBUGGING){

        printf(\
            "%s %s ODE system parameters: \n", \
            identation, function_path \
        );

        printf(\
            "%s \t index: %d \n", \
            identation, (*arg)->index \
        );
        printf(\
            "%s \t equations count: %d \n", \
            identation, (*arg)->eqs_count\
        );
        printf(\
            "%s \t points to save count: %d \n", \
            identation, (*arg)->points_count\
        );
        printf(\
            "%s \t nok init value: %d \n", \
            identation, (*arg)->nok\
        );
        printf(\
            "%s \t nbad init value: %d \n", \
            identation, (*arg)->nbad\
        );
        printf(\
            "%s \t independent variable init value: %.3e\n",\
            identation, (*arg)->x_initial\
        );
        printf(\
            "%s \t independent variable final value: %.3e\n",\
            identation, (*arg)->x_final\
        );

        printf(\
            "%s system initial values \n", identation);

        for(int i=1; i <= (*arg)->eqs_count; i++){
            printf("%s \t y[%d]: %.3e\n", identation, i, (*arg)->y[i]);
        }

        printf(\
            "%s \t name of the system is: %s\n", \
            identation, (*arg)->name_system\
        );
        printf(\
            "%s \t name of the system is: %s\n", \
            identation, (*arg)->name_vars\
        );
    }

    (*arg)->foo = &ode_logistics_foo;

    if(DEBUGGING){
        printf("%s %s ending \n", identation, function_path);
    }

    return;
}

void ode_logistics_free( ODEsystemStruct **arg ){

    const char \
        function_path[] = "All_ODE.c ode_logistics_free", \
        identation[] = "\n\t";

    if(DEBUGGING){
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

void ode_logistics_integrate( ODEsystemStruct **arg ){

    const char function_path[] = "ODE_logistic.c ode_logistics_integrate", \
               identation[] = "\n";

    if(DEBUGGING){
        printf("%s %s starting \n", identation, function_path);
    }

    printf(\
        "%s %s \n\n \
        ODE index [%d] \n \
        \t system name and description - %s \n \
        \t system vars %s \n\n \
        odeint scaling method [%d] \n \
        \t description %s \n\n \
        odeint rkqs method [%d] \n \
        \t description %s \n\n \
        odeint parameters \n \
        \t eps = %.3e \n \
        \t h1 = %.3e \n \
        \t hmin = %.3e \n\n \
        x - start and end point \n \
        \t initial = %.3e \n \
        \t final = %.3e \n\n \
        y - initial conditions \n \
        ", \
        identation, function_path, \
        (*arg)->index, \
        (*arg)->name_system, \
        (*arg)->name_vars, \
        (*arg)->odeint_scaling_method, \
        (*arg)->odeint_scaling_method_description, \
        (*arg)->rkqs_step_method, \
        (*arg)->rkqs_step_method_description, \
        (*arg)->eps, \
        (*arg)->h1, \
        (*arg)->hmin, \
        (*arg)->x_initial, \
        (*arg)->x_final \
    );

    for(int i=1; i <= (*arg)->eqs_count; i++){
        printf("\t\t y[%d] = %.3e \n", i, (*arg)->y[i]);
    }

    double *y;

    y = dvector(1, (*arg)->eqs_count);
    dvector_copy((*arg)->y, y, (*arg)->eqs_count);

    kmax = (*arg)->points_count;
    xp = (*arg)->points_x;
    yp = (*arg)->points_y;

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

#undef DEBUGGING
