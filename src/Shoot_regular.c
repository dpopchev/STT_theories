#include "ExternalHeaders.h"

#define DEBUGGING_shooting_regular_init 0

// how many values of the ODE are given on the left side of the interval
#define KNOWN_LEFT_N 0
#define KNOWN_LEFT_INDEXES 0
#define KNOWN_LEFT_VALUES 0

// how many values of the ODE are unknown and their corresponding guessed values
#define UNKNOWN_LEFT_N 1
#define UNKNOWN_LEFT_INDEXES 1
#define UNKNOWN_LEFT_VALUES 0.8

// how many values of the ODE are given on the right side of the interval
#define KNOWN_RIGHT_N 1
#define KNOWN_RIGHT_INDEXES 1
#define KNOWN_RIGHT_VALUES 1

// how many values of the ODE are unknown and right corresponding guessed values
#define UNKNOWN_RIGHT_N 0
#define UNKNOWN_RIGHT_INDEXES 0
#define UNKNOWN_RIGHT_VALUES 0

void shooting_regular_init(ShootingVarsStruct **arg){

    const char \
        function_path[] = "Shoot_regular.c shooting_regular_init \n", \
        identation[] = "\n";

    if(DEBUGGING_shooting_regular_init){
        printf("%s %s starting \n", identation, function_path);
    }

    (*arg)->known_left_n = KNOWN_LEFT_N;
    if((*arg)->known_left_n){

        int known_left_indexes[] = {
            1e5, // first index is not taken into account
            KNOWN_LEFT_INDEXES
        };

        (*arg)->known_left_indexes = ivector(1,(*arg)->known_left_n);

        double known_left_values[] = {
            1e5, // first index is not taken into account
            KNOWN_LEFT_VALUES
        };

        (*arg)->known_left_values = dvector(1,(*arg)->known_left_n);

        for(int i=1; i<= (*arg)->known_left_n; i++){

            (*arg)->known_left_indexes[i] = known_left_indexes[i];
            (*arg)->known_left_values[i] = known_left_values[i];

            if(DEBUGGING_shooting_regular_init){

                printf(
                    "%s %s known left [%d] = %.3e \n",
                    function_path, identation,
                    (*arg)->known_left_indexes[i], (*arg)->known_left_values[i]
                );
            }

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

    (*arg)->UNknown_left_n = UNKNOWN_LEFT_N;
    if((*arg)->UNknown_left_n){

        int UNknown_left_indexes[] = {
            1e5, // first index is not taken into account
            UNKNOWN_LEFT_INDEXES
        };

        (*arg)->UNknown_left_indexes = ivector(1,(*arg)->UNknown_left_n);

        double UNknown_left_values[] = {
            1e5, // first index is not taken into account
            UNKNOWN_LEFT_VALUES
        };

        (*arg)->UNknown_left_values = dvector(1,(*arg)->UNknown_left_n);

        for(int i=1; i<= (*arg)->known_left_n; i++){

            (*arg)->UNknown_left_indexes[i] = UNknown_left_indexes[i];
            (*arg)->UNknown_left_values[i] = UNknown_left_values[i];

            if(DEBUGGING_shooting_regular_init){

                printf(
                    "%s %s UNknown left [%d] = %.3e \n",
                    function_path, identation,
                    (*arg)->UNknown_left_indexes[i], (*arg)->UNknown_left_values[i]
                );
            }
        }
    }else{
        (*arg)->UNknown_left_indexes = NULL;
        (*arg)->UNknown_left_values = NULL;

        if(DEBUGGING_shooting_regular_init){
            printf(
              "%s %s non UNknown indexes on left side of interval \n",
              function_path, identation
            );
        }

    }

    (*arg)->known_right_n = KNOWN_RIGHT_N;
    if((*arg)->known_right_n){

        int known_right_indexes[] = {
            1e5, // first index is not taken into account
            KNOWN_LEFT_INDEXES
        };

        (*arg)->known_right_indexes = ivector(1,(*arg)->known_right_n);

        double known_right_values[] = {
            1e5, // first index is not taken into account
            KNOWN_LEFT_VALUES
        };

        (*arg)->known_right_values = dvector(1,(*arg)->known_right_n);

        for(int i=1; i<= (*arg)->known_right_n; i++){

            (*arg)->known_right_indexes[i] = known_right_indexes[i];
            (*arg)->known_right_values[i] = known_right_values[i];

            if(DEBUGGING_shooting_regular_init){

                printf(
                    "%s %s known left [%d] = %.3e \n",
                    function_path, identation,
                    (*arg)->known_right_indexes[i], (*arg)->known_right_values[i]
                );
            }

        }
    }else{

        (*arg)->known_right_indexes = NULL;
        (*arg)->known_right_values = NULL;

        if(DEBUGGING_shooting_regular_init){
            printf(
              "%s %s non known indexes on right side of interval \n",
              function_path, identation
            );
        }
    }

    (*arg)->UNknown_right_n = UNKNOWN_RIGHT_N;
    if((*arg)->UNknown_right_n){

        int UNknown_right_indexes[] = {
            1e5, // first index is not taken into account
            UNKNOWN_RIGHT_INDEXES
        };

        (*arg)->UNknown_right_indexes = ivector(1,(*arg)->UNknown_right_n);

        double UNknown_right_values[] = {
            1e5, // first index is not taken into account
            UNKNOWN_RIGHT_VALUES
        };

        (*arg)->UNknown_right_values = dvector(1,(*arg)->UNknown_right_n);

        for(int i=1; i<= (*arg)->UNknown_right_n; i++){

            (*arg)->UNknown_right_indexes[i] = UNknown_right_indexes[i];
            (*arg)->UNknown_right_values[i] = UNknown_right_values[i];

            if(DEBUGGING_shooting_regular_init){

                printf(
                    "%s %s UNknown right [%d] = %.3e \n",
                    function_path, identation,
                    (*arg)->UNknown_right_indexes[i], (*arg)->UNknown_right_values[i]
                );
            }
        }
    }else{

        (*arg)->UNknown_right_indexes = NULL;
        (*arg)->UNknown_right_values = NULL;

        if(DEBUGGING_shooting_regular_init){
            printf(
              "%s %s non UNknown indexes on right side of interval \n",
              function_path, identation
            );
        }
    }

    (*arg)->descriptancy = dvector(1, (*arg)->known_left_n);
    (*arg)->left_achieved = dvector(1, (*arg)->known_left_n);
}

void shooting_regular_free(ShootingVarsStruct **arg){

    free((*arg)->known_left_indexes);
    free((*arg)->UNknown_left_indexes);
    free((*arg)->known_right_indexes);
    free((*arg)->UNknown_right_indexes);
    free((*arg)->known_left_values);
    free((*arg)->UNknown_left_values);
    free((*arg)->known_right_values);
    free((*arg)->UNknown_right_values);
    free((*arg)->descriptancy);
    free((*arg)->left_achieved);
    free((*arg)->right_achieved);
    free((*arg));

    return;
}

void shooting_regular_info_print_stdout(ShootingVarsStruct *arg){

    printf(
        "\n\n Shooting info: \n \
        name: %s \n\n \
        known left side %d: ",
        "shooting regular", arg->known_left_n
    );

    for(int i=1; i <= arg->known_left_n; i++){
        printf("[%d] = %.3e ", arg->known_left_indexes[i], arg->known_left_values[i]);
    }

    printf("\n\t UNknown left side %d: ", arg->UNknown_left_n);

    for(int i=1; i <= arg->UNknown_left_n; i++){
        printf("[%d] = %.3e ", arg->UNknown_left_indexes[i], arg->UNknown_left_values[i]);
    }

    printf("\n\t known right side %d: ", arg->known_right_n);

    for(int i=1; i <= arg->known_right_n; i++){
        printf("[%d] = %.3e ", arg->known_right_indexes[i], arg->known_right_values[i]);
    }

    printf("\n\t UNknown right side %d: ", arg->UNknown_right_n);

    for(int i=1; i <= arg->UNknown_right_n; i++){
        printf("[%d] = %.3e ", arg->UNknown_right_indexes[i], arg->UNknown_right_values[i]);
    }

    return;
}

void shooting_regular_info_print_ResultFile(
  ShootingVarsStruct *arg, FILE *fp
){

    fprintf(
        fp,
        "\n\n Shooting info: \n \
        name: %s \n\n \
        known left side %d: ",
        "shooting regular", arg->known_left_n
    );

    for(int i=1; i <= arg->known_left_n; i++){
        fprintf(fp,"[%d] = %.3e ", arg->known_left_indexes[i], arg->known_left_values[i]);
    }

    fprintf(fp,"\n\t UNknown left side %d: ", arg->UNknown_left_n);

    for(int i=1; i <= arg->UNknown_left_n; i++){
        fprintf(fp,"[%d] = %.3e ", arg->UNknown_left_indexes[i], arg->UNknown_left_values[i]);
    }

    fprintf(fp,"\n\t known right side %d: ", arg->known_right_n);

    for(int i=1; i <= arg->known_right_n; i++){
        fprintf(fp,"[%d] = %.3e ", arg->known_right_indexes[i], arg->known_right_values[i]);
    }

    fprintf(fp,"\n\t UNknown right side %d: ", arg->UNknown_right_n);

    for(int i=1; i <= arg->UNknown_right_n; i++){
        fprintf(fp,"[%d] = %.3e ", arg->UNknown_right_indexes[i], arg->UNknown_right_values[i]);
    }

    return;
}

