#include "ExternalHeaders.h"

#define DEBUGGING_dvector_copy 0

#define PARAM_VAL_CHAR_LENG 16

void dvector_copy( double *copy_from, double *copy_to, int size ){

    const char \
        function_path[] = "General_purpose.c dvector_copy", \
        identation[] = "\n";

    if(DEBUGGING_dvector_copy){
        printf(\
            "%s %s starting to copy dvector with size %d \n", \
            identation, function_path, size \
        );
    }

    for(int i=1; i <= size; i++){

        if(DEBUGGING_dvector_copy){
            printf(\
                "%s %s copy_from[%d] = %.3e copy_to[%d] = %.3e \n", \
                identation, function_path, \
                i, copy_from[i], \
                i, copy_to[i] \
            );
        }

        copy_to[i] = copy_from[i];
    }

    return;
}

FILE* open_file_to_WRITE_ResultFile( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/";
    char **parm_val;

    strcat( full_file_path, arg->name_system );

    parm_val = calloc((size_t)arg->free_parmeters_count_all+1, sizeof(char*));

    for(int i=1; i <= arg->free_parmeters_count_all; i++){

        parm_val[i] = calloc((size_t)PARAM_VAL_CHAR_LENG,sizeof(char));
        sprintf( parm_val[i], "%.3e", arg->free_parmeters_values[i] );
        strcat( full_file_path, "_" );
        strcat( full_file_path, arg->free_parmeters_names[i]);
        strcat( full_file_path, parm_val[i]);

        free(parm_val[i]);
    }

    free(parm_val);

    FILE *fp = fopen(full_file_path,"w");
    if( fp == NULL ){
        printf(
          "\n\n ERROR open_file_to_write_ResultFile oppenning file:\n\n \
          %s ERROR \n\n",
          full_file_path
        );
        perror("Error");

        exit(123);
    }else{
        return fp;
    }

}

FILE* open_file_to_APPEND_ResultFile( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/";
    char **parm_val;

    strcat( full_file_path, arg->name_system );

    parm_val = calloc((size_t)arg->free_parmeters_count_all+1, sizeof(char*));

    for(int i=1; i <= arg->free_parmeters_count_all; i++){

        parm_val[i] = calloc((size_t)PARAM_VAL_CHAR_LENG,sizeof(char));
        sprintf( parm_val[i], "%.3e", arg->free_parmeters_values[i] );
        strcat( full_file_path, "_" );
        strcat( full_file_path, arg->free_parmeters_names[i]);
        strcat( full_file_path, parm_val[i]);

        free(parm_val[i]);
    }

    free(parm_val);

    FILE *fp = fopen(full_file_path,"a");
    if( fp == NULL ){
        printf(
          "\n\n ERROR open_file_to_write_ResultFile oppenning file:\n\n \
          %s ERROR \n\n",
          full_file_path
        );
        perror("Error");

        exit(123);
    }else{
        return fp;
    }

}

FILE* open_file_to_WRITE_LivePlot( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/";
    char **parm_val;

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_live_plot" );

    FILE *fp = fopen(full_file_path,"w");
    if( fp == NULL ){
        printf(
          "\n\n ERROR open_file_to_write_ResultFile oppenning file:\n\n \
          %s ERROR \n\n",
          full_file_path
        );
        perror("Error");

        exit(123);
    }else{
        return fp;
    }

}

FILE* open_file_to_APPEND_LivePlot( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/";
    char **parm_val;

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_live_plot" );

    FILE *fp = fopen(full_file_path,"a");
    if( fp == NULL ){
        printf(
          "\n\n ERROR open_file_to_write_ResultFile oppenning file:\n\n \
          %s ERROR \n\n",
          full_file_path
        );
        perror("Error");

        exit(123);
    }else{
        return fp;
    }

}

#undef PARAM_VAL_CHAR_LENG
#undef DEBUGGING_dvector_copy
