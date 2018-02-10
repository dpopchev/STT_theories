#include "ExternalHeaders.h"

#define DEBUGGING_dvector_copy 0
#define DEBUGGING_dvector_copy_to_index 0
#define DEBUGGING_OPEN_WRITE_RESULT 0

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

void dvector_copy_to_index( double *copy_from, double *copy_to, int size, int *indexes){

    const char \
        function_path[] = "General_purpose.c dvector_copy_to_index", \
        identation[] = "\n";

    if(DEBUGGING_dvector_copy_to_index){
        printf(\
            "%s %s starting to copy dvector with size %d \n", \
            identation, function_path, size \
        );
    }

    for(int i=1; i <= size; i++){

        if(DEBUGGING_dvector_copy){
            printf(
                "%s %s copy_from[%d] = %.3e copy_to[%d] = %.3e \n",
                identation, function_path,
                indexes[i], copy_from[indexes[i]],
                indexes[i], copy_to[indexes[i]]
            );
        }

        copy_to[indexes[i]] = copy_from[indexes[i]];
    }

    return;
}

FILE* open_file_to_WRITE_ResultFile( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/";
    char **parm_val;

    if(DEBUGGING_OPEN_WRITE_RESULT){
        printf("\n General_purpose.c open_file_to_WRITE_ResultFile \n");
    }

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_");
    strcat( full_file_path, arg->eoseq_name);

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
          %s ERROR, line 104\n\n",
          full_file_path
        );
        perror("Error");

        exit(104);
    }else{
        return fp;
    }
}

FILE* open_file_to_APPEND_ResultFile( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/";
    char **parm_val;

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_");
    strcat( full_file_path, arg->eoseq_name);

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
          %s ERROR line 141 \n\n",
          full_file_path
        );
        perror("Error");

        exit(141);
    }else{
        return fp;
    }
}

FILE* open_file_to_WRITE_LivePlot( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/live_plot_";
    char **parm_val;

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_");
    strcat( full_file_path, arg->eoseq_name);

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
          %s ERROR line 178\n\n",
          full_file_path
        );
        perror("Error");

        exit(178);
    }else{
        return fp;
    }
}

FILE* open_file_to_APPEND_LivePlot( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/live_plot_";
    char **parm_val;

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_");
    strcat( full_file_path, arg->eoseq_name);

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
          %s ERROR line 215 \n\n",
          full_file_path
        );
        perror("Error");

        exit(215);
    }else{
        return fp;
    }
}

FILE* open_file_to_WRITE_LivePlot_solver( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/live_plot_solver_";
    char **parm_val;

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_");
    strcat( full_file_path, arg->eoseq_name);

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
          %s ERROR line 178\n\n",
          full_file_path
        );
        perror("Error");

        exit(178);
    }else{
        return fp;
    }
}

FILE* open_file_to_APPEND_LivePlot_solver( ODEsystemStruct *arg ){

    char full_file_path[128]="/home/dimitar/projects/STT_theories/results/live_plot_solver_";
    char **parm_val;

    strcat( full_file_path, arg->name_system );
    strcat( full_file_path, "_");
    strcat( full_file_path, arg->eoseq_name);

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
          %s ERROR line 215 \n\n",
          full_file_path
        );
        perror("Error");

        exit(215);
    }else{
        return fp;
    }
}

#undef PARAM_VAL_CHAR_LENG
#undef DEBUGGING_dvector_copy
#undef DEBUGGING_OPEN_WRITE_RESULT
#undef DEBUGGING_dvector_copy_to_index
