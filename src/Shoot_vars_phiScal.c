#include "ExternalHeaders.h"

#define DEBUGGING_shooting_phiScal_regular_init 0

// how many values of the ODE are given on the left side of the interval
#define KNOWN_LEFT_N 2
#define KNOWN_LEFT_INDEXES 2,4
#define KNOWN_LEFT_VALUES 0,0

// how many values of the ODE are unknown and their corresponding guessed values
#define UNKNOWN_LEFT_N 1
#define UNKNOWN_LEFT_INDEXES 1
#define UNKNOWN_LEFT_VALUES -5e-2

// how many values of the ODE are given on the right side of the interval
#define KNOWN_RIGHT_N 1
#define KNOWN_RIGHT_INDEXES 1
#define KNOWN_RIGHT_VALUES 0

// how many values of the ODE are unknown and right corresponding guessed values
#define UNKNOWN_RIGHT_N 2
#define UNKNOWN_RIGHT_INDEXES 2,4
#define UNKNOWN_RIGHT_VALUES 0,0

// TODO this is not ready for the fitting point
#define SHOOT_FREE_N 2
#define SHOOT_FREE_INDEXES 3,5

void shooting_phiScal_init(ShootingVarsStruct **arg){

    const char \
        function_path[] = "Shoot_vars.c shooting_phiScal_regular_init ", \
        identation[] = "\n\t";

    if(DEBUGGING_shooting_phiScal_regular_init){
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

            if(DEBUGGING_shooting_phiScal_regular_init){

                printf(
                    "%s %s known left [%d] = %.3e \n",
                    identation, function_path,
                    (*arg)->known_left_indexes[i], (*arg)->known_left_values[i]
                );
            }
        }
    }else{

        (*arg)->known_left_indexes = NULL;
        (*arg)->known_left_values = NULL;

        if(DEBUGGING_shooting_phiScal_regular_init){
            printf(
              "%s %s non known indexes on left side of interval \n",
              identation, function_path
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

        for(int i=1; i<= (*arg)->UNknown_left_n; i++){

            (*arg)->UNknown_left_indexes[i] = UNknown_left_indexes[i];
            (*arg)->UNknown_left_values[i] = UNknown_left_values[i];

            if(DEBUGGING_shooting_phiScal_regular_init){

                printf(
                    "%s %s (init guess) UNknown left [%d] = %.3e \n",
                    identation, function_path,
                    (*arg)->UNknown_left_indexes[i], (*arg)->UNknown_left_values[i]
                );
            }
        }
    }else{
        (*arg)->UNknown_left_indexes = NULL;
        (*arg)->UNknown_left_values = NULL;

        if(DEBUGGING_shooting_phiScal_regular_init){
            printf(
              "%s %s non UNknown indexes on left side of interval \n",
              identation, function_path
            );
        }
    }

    (*arg)->known_right_n = KNOWN_RIGHT_N;
    if((*arg)->known_right_n){

        int known_right_indexes[] = {
            1e5, // first index is not taken into account
            KNOWN_RIGHT_INDEXES
        };

        (*arg)->known_right_indexes = ivector(1,(*arg)->known_right_n);

        double known_right_values[] = {
            1e5, // first index is not taken into account
            KNOWN_RIGHT_VALUES
        };

        (*arg)->known_right_values = dvector(1,(*arg)->known_right_n);

        for(int i=1; i<= (*arg)->known_right_n; i++){

            (*arg)->known_right_indexes[i] = known_right_indexes[i];
            (*arg)->known_right_values[i] = known_right_values[i];

            if(DEBUGGING_shooting_phiScal_regular_init){

                printf(
                    "%s %s known right [%d] = %.3e \n",
                    identation, function_path,
                    (*arg)->known_right_indexes[i], (*arg)->known_right_values[i]
                );
            }
        }
    }else{

        (*arg)->known_right_indexes = NULL;
        (*arg)->known_right_values = NULL;

        if(DEBUGGING_shooting_phiScal_regular_init){
            printf(
              "%s %s non known indexes on right side of interval \n",
              identation, function_path
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

            if(DEBUGGING_shooting_phiScal_regular_init){

                printf(
                    "%s %s (init guesses) UNknown right [%d] = %.3e \n",
                    identation, function_path,
                    (*arg)->UNknown_right_indexes[i], (*arg)->UNknown_right_values[i]
                );
            }
        }
    }else{

        (*arg)->UNknown_right_indexes = NULL;
        (*arg)->UNknown_right_values = NULL;

        if(DEBUGGING_shooting_phiScal_regular_init){
            printf(
              "%s %s non (init guesses) UNknown indexes on right side of interval \n",
              identation, function_path
            );
        }
    }

    (*arg)->newt_n = (*arg)->UNknown_left_n + (*arg)->UNknown_right_n;
    (*arg)->newt_v = dvector(1, (*arg)->newt_n);
    (*arg)->newt_f = dvector(1, (*arg)->newt_n);

    if((*arg)->newt_n){

        for(int i=1; i<= (*arg)->UNknown_left_n; i++){

            (*arg)->newt_v[i] = (*arg)->UNknown_left_values[i];

            if(DEBUGGING_shooting_phiScal_regular_init){

                printf(
                    "%s %s newt_v[%d] = %.3e \t newt_f[%d] = %.3e \n",
                    identation, function_path,
                    i, (*arg)->newt_v[i],
                    i, (*arg)->newt_f[i]
                );
            }
        }

        for(int i = 1; i<= (*arg)->UNknown_right_n; i++){

            (*arg)->newt_v[(*arg)->UNknown_left_n + i] = (*arg)->UNknown_right_values[i];

            if(DEBUGGING_shooting_phiScal_regular_init){

                printf(
                    "%s %s newt_v[%d] = %.3e \t newt_f[%d] = %.3e \n",
                    identation, function_path,
                    (*arg)->UNknown_left_n + i, (*arg)->newt_v[(*arg)->UNknown_left_n + i],
                    (*arg)->UNknown_left_n + i, (*arg)->newt_f[(*arg)->UNknown_left_n + i]
                );
            }
        }
    }else{
        printf(
            "%s %s Here we should have something different, line 237 \n",
            identation, function_path
        );

        exit(237);
    }

    (*arg)->descriptancy = dvector(1, (*arg)->newt_n);
    (*arg)->left_achieved = dvector(1, (*arg)->UNknown_left_n);
    (*arg)->right_achieved = dvector(1, (*arg)->UNknown_right_n);

    (*arg)->shoot_free_n = SHOOT_FREE_N;
    if((*arg)->shoot_free_n){

        int shoot_free_indexes[] = {
            1e5, // first index is not taken into account
            SHOOT_FREE_INDEXES
        };

        (*arg)->shoot_free_indexes = ivector(1, (*arg)->shoot_free_n);

        for(int i=1; i <=(*arg)->shoot_free_n; i++){

            (*arg)->shoot_free_indexes[i] = shoot_free_indexes[i];

            if(DEBUGGING_shooting_phiScal_regular_init){

                printf(
                    "%s %s shoot free index [%d] \n",
                    identation, function_path,
                    (*arg)->shoot_free_indexes[i]
                );
            }
        }
    }

    return;
}

void shooting_phiScal_free(ShootingVarsStruct **arg){

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
    free((*arg)->newt_v);
    free((*arg)->newt_f);
    free((*arg)->shoot_free_indexes);
    free((*arg));

    return;
}

void shooting_phiScal_info_print_stdout(ShootingVarsStruct *arg){

    printf(
        "\n shooting_phiScal info: \n"
        "\n\t General info: \n"
        "\t\t name %s \n",
        "general shooting_phiScal, pay attention"
    );

    printf(
        "\n\t Known left side\n"
        "\t\t count %d -> ",
        arg->known_left_n
    );

    for(int i=1; i <= arg->known_left_n; i++){
        printf(
            "\n\t\t\t [%d] = %.3e ",
            arg->known_left_indexes[i], arg->known_left_values[i]
        );
    }

    printf(
        "\n\n\t UNknown left side (init guess) \n"
        "\t\t count %d -> ",
        arg->UNknown_left_n
    );

    for(int i=1; i <= arg->UNknown_left_n; i++){
        printf(
            "\n\t\t\t [%d] = %.3e ",
            arg->UNknown_left_indexes[i], arg->UNknown_left_values[i]
        );
    }

    printf(
        "\n\n\t Known right side \n"
        "\t\t count %d -> ",
        arg->known_right_n
    );

    for(int i=1; i <= arg->known_right_n; i++){
        printf(
            "\n\t\t\t [%d] = %.3e ",
            arg->known_right_indexes[i], arg->known_right_values[i]
        );
    }

    printf(
        "\n\n\t UNknown right side (init guess) \n"
        "\t\t count %d -> ",
        arg->UNknown_right_n
    );

    for(int i=1; i <= arg->UNknown_right_n; i++){
        printf(
            "\n\t\t\t [%d] = %.3e ",
            arg->UNknown_right_indexes[i], arg->UNknown_right_values[i]
        );
    }

    printf(
        "\n\n\t Shoot free indexes: \n"
        "\t\t count %d -> ",
        arg->shoot_free_n
    );

    for(int i=1; i <= arg->shoot_free_n; i++){
        printf(
            "\n\t\t\t [%d]",
            arg->shoot_free_indexes[i]
        );
    }

    printf("\n");
    return;
}

void shooting_phiScal_info_print_ResultFile( ShootingVarsStruct *arg, FILE *fp ){

    fprintf(
        fp,
        "\n shooting_phiScal info: \n"
        "\n\t General info: \n"
        "\t\t name %s \n",
        "general shooting_phiScal, pay attention"
    );

    fprintf(
        fp,
        "\n\t Known left side\n"
        "\t\t count %d -> ",
        arg->known_left_n
    );

    for(int i=1; i <= arg->known_left_n; i++){
        fprintf(
            fp,
            "\n\t\t\t [%d] = %.3e ",
            arg->known_left_indexes[i], arg->known_left_values[i]
        );
    }

    fprintf(
        fp,
        "\n\n\t UNknown left side (init guess) \n"
        "\t\t count %d -> ",
        arg->UNknown_left_n
    );

    for(int i=1; i <= arg->UNknown_left_n; i++){
        fprintf(
            fp,
            "\n\t\t\t [%d] = %.3e ",
            arg->UNknown_left_indexes[i], arg->UNknown_left_values[i]
        );
    }

    fprintf(
        fp,
        "\n\n\t known right side \n"
        "\t\t count %d -> ",
        arg->known_right_n
    );

    for(int i=1; i <= arg->known_right_n; i++){
        fprintf(
            fp,
            "\n\t\t\t [%d] = %.3e ",
            arg->known_right_indexes[i], arg->known_right_values[i]
        );
    }

    fprintf(
        fp,
        "\n\n\t UNknown right side (init guess) \n"
        "\t\t count %d -> ",
        arg->UNknown_right_n
    );

    for(int i=1; i <= arg->UNknown_right_n; i++){
        fprintf(
            fp,
            "\n\t\t\t [%d] = %.3e ",
            arg->UNknown_right_indexes[i], arg->UNknown_right_values[i]
        );
    }

    fprintf(fp,"\n\n");

    //TODO this better not be here
    fprintf( fp, "#p_c phiScal_c M AR rho_c delta_phiScal\n");

    fclose(fp);
    return;
}

void shooting_phiScal_check(ShootingVarsStruct *arg_shoot, ODEsystemStruct *arg_ode){

    if(
      arg_shoot->known_left_n + arg_shoot->UNknown_left_n != arg_ode->eqs_count - SHOOT_FREE_N
    ){
        printf(
          "\n Shoot_vars.c %d + %d != %d line 344 \n",
          arg_shoot->known_left_n, arg_shoot->UNknown_left_n,
          arg_ode->eqs_count
        );

        exit(344);
    }else if(
      arg_shoot->known_right_n + arg_shoot->UNknown_right_n != arg_ode->eqs_count - SHOOT_FREE_N
    ){
        printf(
          "\n Shoot_vars.c %d + %d != %d line 354 \n",
          arg_shoot->known_right_n, arg_shoot->UNknown_right_n,
          arg_ode->eqs_count
        );

        exit(354);
    }else if(
      arg_shoot->known_left_n + arg_shoot->known_right_n != arg_ode->eqs_count - SHOOT_FREE_N
    ){
        printf(
          "\n Shoot_vars.c %d + %d != %d line 364 \n",
          arg_shoot->known_left_n, arg_shoot->known_right_n,
          arg_ode->eqs_count
        );

        exit(364);
    }else if(
      arg_shoot->UNknown_left_n + arg_shoot->UNknown_right_n != arg_ode->eqs_count - SHOOT_FREE_N
    ){
        printf(
          "\n Shoot_vars.c %d + %d != %d line 374 \n",
          arg_shoot->UNknown_left_n, arg_shoot->UNknown_right_n,
          arg_ode->eqs_count
        );

        exit(374);

    }else if( arg_shoot->known_left_n != arg_shoot->UNknown_right_n ){
        printf(
          "\n Shoot_vars.c known_left %d != %d UNknown_right line 382 \n",
          arg_shoot->known_left_n, arg_shoot->UNknown_right_n
        );

        exit(382);
    }else if( arg_shoot->UNknown_left_n != arg_shoot->known_right_n ){
        printf(
          "\n Shoot_vars.c UNknown_left %d != %d known_right line 389 \n",
          arg_shoot->UNknown_left_n, arg_shoot->known_right_n
        );

        exit(389);
    }

    return;
}

#undef DEBUGGING_shooting_phiScal_regular_init
#undef KNOWN_LEFT_N
#undef KNOWN_LEFT_INDEXES
#undef KNOWN_LEFT_VALUES
#undef UNKNOWN_LEFT_N
#undef UNKNOWN_LEFT_INDEXES
#undef UNKNOWN_LEFT_VALUES
#undef KNOWN_RIGHT_N
#undef KNOWN_RIGHT_INDEXES
#undef KNOWN_RIGHT_VALUES
#undef UNKNOWN_RIGHT_N
#undef UNKNOWN_RIGHT_INDEXES
#undef UNKNOWN_RIGHT_VALUES
