#include <ExternalHeaders.h>

#define DEBUGGING 2
#define DEBUGGING_units 2
#define DEBUGGING_epsilon 2

// sets the units factors as we are working dimensionless
static void calculate_gv_units( void ){

    const char function_path[] = "main.c calculate_gv_units", \
               place_holder[] = "\n\t";

    if(DEBUGGING_units){
        printf("%s %s starting \n", place_holder, function_path);
    }

    GV_RADCORD_UNITS = 1e-3*GV_G*GV_MSUN/pow(GV_C,2);
    if(DEBUGGING_units >= 2){
        printf("%s %s GV_RADCORD_UNITS = %.2e \n", place_holder, function_path, GV_RADCORD_UNITS);
    }

    GV_DENSITY_UNITS = 1e-3*pow(GV_C,6)/(pow(GV_G,3)*pow(GV_MSUN,2));
    if(DEBUGGING_units >= 2){
        printf("%s %s GV_DENSITY_UNITS = %.2e \n", place_holder, function_path, GV_DENSITY_UNITS);
    }

    GV_PRESSURE_UNITS = pow(GV_C,8)/(pow(GV_G,3)*pow(GV_MSUN,2))*10.0;
    if(DEBUGGING_units >= 2){
        printf("%s %s GV_PRESSURE_UNITS = %.2e \n", place_holder, function_path, GV_PRESSURE_UNITS);
    }

    if(DEBUGGING_units){
        printf("%s %s ending \n", place_holder, function_path);
    }

    return;
}

static void get_machine_eps( void ){

    const char function_path[] = "main.c get_machine_eps", \
               place_holder[] = "\n\t";

    if(DEBUGGING_epsilon){
        printf("%s %s starting \n", place_holder, function_path);
    }

    GV_MACHINE_EPSILON = 1;
    while( 1 + GV_MACHINE_EPSILON/2 != 1 ){
        GV_MACHINE_EPSILON /= 2;
    }

    GV_MACHINE_EPSILON_CUBROOT = cbrt(GV_MACHINE_EPSILON);
    GV_MACHINE_EPSILON_SQRT = sqrt(GV_MACHINE_EPSILON);

    if(DEBUGGING_epsilon){
        printf(
            "%s %s \n\t\t eps = %.2e \n\t\t cbrt(eps) = %.2e \n\t\t sqrt(eps) = %.2e  \n", \
            place_holder, function_path, GV_MACHINE_EPSILON, \
            GV_MACHINE_EPSILON_CUBROOT, GV_MACHINE_EPSILON_SQRT \
        );
    }

    if(DEBUGGING_epsilon){
        printf("%s %s ending \n", place_holder, function_path);
    }

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

    return 0;
}
