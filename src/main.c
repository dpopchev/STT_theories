#include "ExternalHeaders.h"

#define DEBUGGING 1

#define ODE_SYSTEM_INDEX 1
#define ODE_SYSTEM_EQS_COUNT 3
#define ODE_SYSTEM_POINTS_COUNT 0

// the constants of nature
double \
    // mass of the sun in kg
    GV_MSUN = 1.9891e30, \
    // gravitational constant in  m^3*kg^-1*s^-2
    GV_G = 6.67384e-11, \
    // speed of light in m*s^-1
    GV_C = 299792458, \
    // FIXME: What is GV_MBAR constant?
    GV_MBAR = 1.66e-24, \
    // Nice pi number
    GV_PI = 3.14159265358979323846, \
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

// name and simple description of the ODE system
char *ODEsystem_names[] = { \
        "System 1, simpliest", \
        "System 2, not so simple " \
    };
char *ODEsystem_vars[] = { \
        "phi, etc..", \
        "phi and other variables" \
    };

// sets the units factors as we are working dimensionless
static void calculate_gv_units( void ){

    const char function_path[] = "main.c calculate_gv_units", \
               place_holder[] = "\n\t";

    if(DEBUGGING){
        printf("%s %s starting \n", place_holder, function_path);
    }

    GV_RADCORD_UNITS = 1e-3*GV_G*GV_MSUN/pow(GV_C,2);
    if(DEBUGGING >= 2){
        printf("%s %s GV_RADCORD_UNITS = %.2e \n", place_holder, function_path, GV_RADCORD_UNITS);
    }

    GV_DENSITY_UNITS = 1e-3*pow(GV_C,6)/(pow(GV_G,3)*pow(GV_MSUN,2));
    if(DEBUGGING >= 2){
        printf("%s %s GV_DENSITY_UNITS = %.2e \n", place_holder, function_path, GV_DENSITY_UNITS);
    }

    GV_PRESSURE_UNITS = pow(GV_C,8)/(pow(GV_G,3)*pow(GV_MSUN,2))*10.0;
    if(DEBUGGING >= 2){
        printf("%s %s GV_PRESSURE_UNITS = %.2e \n", place_holder, function_path, GV_PRESSURE_UNITS);
    }

    if(DEBUGGING){
        printf("%s %s ending \n", place_holder, function_path);
    }

    return;
}

static void get_machine_eps( void ){

    const char function_path[] = "main.c get_machine_eps", \
               place_holder[] = "\n\t";

    if(DEBUGGING){
        printf("%s %s starting \n", place_holder, function_path);
    }

    GV_MACHINE_EPSILON = 1;
    while( 1 + GV_MACHINE_EPSILON/2 != 1 ){
        GV_MACHINE_EPSILON /= 2;
    }

    GV_MACHINE_EPSILON_CUBROOT = cbrt(GV_MACHINE_EPSILON);
    GV_MACHINE_EPSILON_SQRT = sqrt(GV_MACHINE_EPSILON);

    if(DEBUGGING){
        printf(
            "%s %s \n\t\t eps = %.2e \n\t\t cbrt(eps) = %.2e \n\t\t sqrt(eps) = %.2e  \n", \
            place_holder, function_path, GV_MACHINE_EPSILON, \
            GV_MACHINE_EPSILON_CUBROOT, GV_MACHINE_EPSILON_SQRT \
        );
    }

    if(DEBUGGING){
        printf("%s %s ending \n", place_holder, function_path);
    }

    return;
}

static void ODEsystem_init( ODEsystemStruct **self ){

    const char function_path[] = "main.c ODEsystem_init", \
               place_holder[] = "\n\t";

    if(DEBUGGING){
        printf("%s %s starting \n", place_holder, function_path);
    }

    if( ( *self = calloc(1,sizeof(ODEsystemStruct)) ) == NULL ){
        printf("%s %s error in calloc for ode \n", place_holder, function_path);
        exit(123);
    };

    (*self)->index = ODE_SYSTEM_INDEX;
    (*self)->eqs_count = ODE_SYSTEM_EQS_COUNT;
    (*self)->points_count = ODE_SYSTEM_POINTS_COUNT;
    (*self)->nok = 0;
    (*self)->nbad = 0;
    (*self)->x = 0;
    (*self)->y = dvector(1, (*self)->eqs_count);

    if((*self)->points_count){
        (*self)->points_x = dvector(1, (*self)->points_count);
        (*self)->points_y = dmatrix(1, (*self)->eqs_count, 1, (*self)->points_count);
    }

    strcpy((*self)->name_system, ODEsystem_names[(*self)->index]);
    strcpy((*self)->name_vars, ODEsystem_vars[(*self)->index]);

    if(DEBUGGING){

        printf("%s %s ODE system parameters: \n", place_holder, function_path);

        printf("\t\t index: %d \n", (*self)->index);
        printf("\t\t equations count: %d \n", (*self)->eqs_count);
        printf("\t\t points to save count: %d \n", (*self)->points_count);
        printf("\t\t nok init value: %d \n", (*self)->nok);
        printf("\t\t nbad init value: %d \n", (*self)->nbad);
        printf("\t\t independent variable init value: %.3e\n", (*self)->x);

        for(int i=1; i <= (*self)->eqs_count; i++){
            printf("\t\t y[%d]: %.3e\n", i, (*self)->y[i]);
        }

        printf("\t\t name of the system is: %s\n", (*self)->name_system);
        printf("\t\t name of the system is: %s\n", (*self)->name_vars);

    }

    if(DEBUGGING){
        printf("%s %s ending \n", place_holder, function_path);
    }

    return;
}

static void free_ODEsystem( ODEsystemStruct **self ){

    const char function_path[] = "main.c free_ODEsystem", \
               place_holder[] = "\n\t";

    if(DEBUGGING){
        printf("%s %s starting \n", place_holder, function_path);
    }

    if( (*self)->eqs_count ){
        free_dvector( (*self)->y, 1, (*self)->eqs_count );
    }else{
        free( (*self)->y );
    }

    if( (*self)->points_count ){
        free_dvector( (*self)->points_x, 1, (*self)->points_count );
        free_dmatrix( \
            (*self)->points_y, \
            1, (*self)->eqs_count, \
            1, (*self)->points_count \
        );
    }else{
        free((*self)->points_x);
        free((*self)->points_y);
    }

    free(*self);

    return;
}

int main( void ){

    const char function_path[] = "main.c get_machine_eps", \
               place_holder[] = "\n";

    if(DEBUGGING){
        printf("%s %s starting \n", place_holder, function_path);
    }

    calculate_gv_units();
    get_machine_eps();

    ODEsystemStruct *ode;

    ODEsystem_init( &ode );

    if(DEBUGGING){
        printf("%s %s ending \n", place_holder, function_path);
    }

    free_ODEsystem(&ode);

    return 0;
}
