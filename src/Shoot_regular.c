#include "ExternalHeaders.h"

#define DEBUGGING_shooting_regular_init

// how many values of the ODE are given on the left side of the interval
#define KNOWN_LEFT_N 0
#define KNOWN_LEFT_INDEXES 0
#define KNOWN_LEFT_VALUES 0

void shooting_regular_init(ShootingVarsStruct **arg){

    const char \
        function_path[] = "Shoot_regular.c shooting_regular_init \n", \
        identation[] = "\n";

    if(DEBUGGING_shooting_regular_init){
        printf("%s %s starting \n", function_path, identation);
    }

    (*arg)->known_left_n = KNOWN_LEFT_N;
    if((*arg)->known_left_n){

        int known_left_indexes[] = {
            1e5, // first index is not taken into account
            KNOWN_LEFT_INDEXES
        }

        double known_left_values[] = {
            1e5, // first index is not taken into account
            KNOWN_LEFT_VALUES
        }

        (*arg)->known_left_indexes = ivector(1,(*arg)->known_left_n);
        (*arg)->known_left_values = dvector(1,(*arg)->known_left_n);

        for(int i=1; i<= (*arg)->known_left_n; i++){

        }

    }else{

        (*arg)->known_left_indexes = NULL;
        (*arg)->known_left_values = NULL;

        if(DEBUGGING_shooting_regular_init){
            printf(
              "%s %s non known indexes on left side of interval \n",
              function_path, identation
            );
        }
    }


}
