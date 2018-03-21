#include "ExternalHeaders.h"

#define DEBUGGING 0

// 0 - ode logistics, for test purposes
#define ODE_SYSTEM_INDEX 1

// the constants of nature
double \
    // mass of the sun in kg
    GV_MSUN = 1.9891e30, \
    // gravitational constant in  m^3*kg^-1*s^-2
    GV_G = 6.67384e-11, \
    // speed of light in m*s^-1
    GV_C = 299792458, \
    // Barionic mass for ESOII to calculate something of the sun about bounding energy
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
    GV_MACHINE_EPSILON_CUBROOT, \
    // units for moment of inertia gcm^2
    GV_J_UNITS;

// sets the units factors as we are working dimensionless
static void calculate_gv_units( void ){

    const char function_path[] = "main.c calculate_gv_units", \
               identation[] = "\n\t";

    if(DEBUGGING){
        printf("%s %s starting \n", identation, function_path);
    }

    GV_RADCORD_UNITS = 1e-3*GV_G*GV_MSUN/pow(GV_C,2);
    if(DEBUGGING >= 2){
        printf(\
            "%s %s GV_RADCORD_UNITS = %.2e \n", \
            identation, function_path, GV_RADCORD_UNITS\
        );
    }

    GV_DENSITY_UNITS = 1e-3*pow(GV_C,6)/(pow(GV_G,3)*pow(GV_MSUN,2));
    if(DEBUGGING >= 2){
        printf(\
            "%s %s GV_DENSITY_UNITS = %.2e \n", \
            identation, function_path, GV_DENSITY_UNITS\
        );
    }

    GV_PRESSURE_UNITS = pow(GV_C,8)/(pow(GV_G,3)*pow(GV_MSUN,2))*10.0;
    if(DEBUGGING >= 2){
        printf(\
            "%s %s GV_PRESSURE_UNITS = %.2e \n", \
            identation, function_path, GV_PRESSURE_UNITS \
        );
    }

    GV_J_UNITS = 1e-1*pow(GV_G,2)*pow(GV_MSUN,3)/pow(GV_C,4);

    if(DEBUGGING){
        printf("%s %s ending \n", identation, function_path);
    }

    return;
}

static void get_machine_eps( void ){

    const char function_path[] = "main.c get_machine_eps", \
               identation[] = "\n\t";

    if(DEBUGGING){
        printf("%s %s starting \n", identation, function_path);
    }

    GV_MACHINE_EPSILON = 1;
    while( 1 + GV_MACHINE_EPSILON/2 != 1 ){
        GV_MACHINE_EPSILON /= 2;
    }

    GV_MACHINE_EPSILON_CUBROOT = cbrt(GV_MACHINE_EPSILON);
    GV_MACHINE_EPSILON_SQRT = sqrt(GV_MACHINE_EPSILON);

    if(DEBUGGING){
        printf(
            "%s %s \n\t\t eps = %.2e \n\t\t cbrt(eps) = %.2e \n\t\t sqrt(eps) = %.2e \n", \
            identation, function_path, GV_MACHINE_EPSILON, \
            GV_MACHINE_EPSILON_CUBROOT, GV_MACHINE_EPSILON_SQRT \
        );
    }

    if(DEBUGGING){
        printf("%s %s ending \n", identation, function_path);
    }

    return;
}

int main( void ){

    const char \
        function_path[] = "main.c", \
        identation[] = "\n",
        *ODE_equtions_list[] = {
            "logistic equation"
        };

    if(DEBUGGING){
        printf("%s %s starting \n", identation, function_path);
    }

    srand((unsigned int)time(NULL));

    int ODE_index = ODE_SYSTEM_INDEX;

    calculate_gv_units();
    get_machine_eps();

    // single integration of the phiScal
    //single_integrate_phiScal();

    // single shooting of for the phiScal
    //single_shoot_regular();

    // iterate through several central pressures
    //iterate_pressure_shoot_regular();

    //single shoot and iterate the infinity
    //single_shoot_regular_iterate_inf();

    // iterate through central pressure and infinity
    //iterate_pressure_inf_shoot_regular();

    // start using the phiScal to calculate the Omega
    //single_shoot_regular_phiScal_J();

    // start using the phiScal to calculate the Omega
    // as above but with iteration for the inf
    //single_shoot_regular_phiScal_J_iterate_inf();

    // start using the phiScal to calculate the Omega
    // as above but with iteration for the inf
    // iterating over pressures
    single_shoot_regular_phiScal_J_iterate_inf_iterpres();

    //ODEsystemStruct *ode;

    //switch(ODE_index){

        //case 1:
            //if(DEBUGGING){
                //printf(
                    //"%s %s ODE_index = %d to init, %s \n",
                    //identation,
                    //function_path,
                    //ODE_index,
                    //ODE_equtions_list[ODE_index]
                //);
            //}

            //ode_phiScal_init( &ode );

            //ode_phiScal_compute_parameters(ode);

            //break;

        //case 2:
            //if(DEBUGGING){
                //printf(
                    //"%s %s ODE_index = %d to init, %s \n",
                    //identation,
                    //function_path,
                    //ODE_index,
                    //ODE_equtions_list[ODE_index]
                //);
            //}

            //ode_phiScal_I_init( &ode );

            //ode_phiScal_I_compute_parameters(ode);

            //break;

        //default:
            //printf(
                //"%s %s ODE_index = %d UNKNOWN, exiting.., line 153 \n",
                //identation,
                //function_path,
                //ODE_index
            //);
            //exit(153);
    //}

    //if(DEBUGGING){
        //printf("%s %s ending \n", identation, function_path);
    //}

    //switch(ODE_index){

        //case 1 :
            //if(DEBUGGING){
                //printf(
                    //"%s %s ODE_index = %d to free, %s \n",
                    //identation,
                    //function_path,
                    //ODE_index,
                    //ODE_equtions_list[ODE_index]
                //);
            //}

            //ode_phiScal_free( &ode );
            //break;

        //case 2 :
            //if(DEBUGGING){
                //printf(
                    //"%s %s ODE_index = %d to free, %s \n",
                    //identation,
                    //function_path,
                    //ODE_index,
                    //ODE_equtions_list[ODE_index]
                //);
            //}

            //ode_phiScal_I_free( &ode );
            //break;

        //default :
            //printf(
                //"%s %s ODE_index = %d UNKNOWN, exiting.. , line 183\n",
                //identation,
                //function_path,
                //ODE_index
            //);
            //exit(183);
    //}

    return 0;
}

#undef ODE_SYSTEM_INDEX
#undef DEBUGGING
