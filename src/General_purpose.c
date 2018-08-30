#include "ExternalHeaders.h"

#define DEBUGGING_dvector_copy 0
#define DEBUGGING_dvector_copy_to_index 0
#define DEBUGGING_dvector_copy_to_index_opps 0
#define DEBUGGING_OPEN_WRITE_RESULT 0

#define PARAM_VAL_CHAR_LENG 16

void dvector_copy( double *src, double *dst, int size ){

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
                "%s %s src[%d] = %e will go to dst[%d] = %e \n", \
                identation, function_path, \
                i, src[i], \
                i, dst[i] \
            );
        }

        dst[i] = src[i];
    }

    return;
}

void dvector_copy_to_index( double *src, double *dst, int size, int *indexes){

    const char
        function_path[] = "General_purpose.c dvector_copy_to_index",
        identation[] = "\n";

    if(DEBUGGING_dvector_copy_to_index){
        printf(
            "%s %s starting to copy dvector with size %d \n",
            identation, function_path, size
        );
    }

    for(int i=1; i <= size; i++){

        if(DEBUGGING_dvector_copy_to_index){
            printf(
                "%s %s src[%d] = %e will go to dst[%d] = %e \n",
                identation, function_path,
                i, src[i],
                indexes[i], dst[indexes[i]]
            );
        }

        dst[indexes[i]] = src[i];
    }

    return;
}

void dvector_copy_to_index_opps( double *src, double *dst, int size, int *indexes){

    const char \
        function_path[] = "General_purpose.c dvector_copy_to_index_opps",
        identation[] = "\n";

    if(DEBUGGING_dvector_copy_to_index_opps){
        printf(\
            "%s %s starting to copy dvector with size %d \n",
            identation, function_path, size
        );
    }

    for(int i=1; i <= size; i++){

        if(DEBUGGING_dvector_copy_to_index_opps){
            printf(
                "%s %s src[%d] = %.3e will go to dst[%d] = %.3e \n",
                identation, function_path,
                indexes[i], src[indexes[i]],
                i, dst[i]
            );
        }

        dst[i] = src[indexes[i]];
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

    //printf("\n %s \n", full_file_path);

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

    //printf("\n %s \n", full_file_path);

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

    //printf("\n %s \n", full_file_path);

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

    //printf("\n %s \n", full_file_path);

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

void LivePlot_open(const char *ode_name, char *eos_name, double *pars){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/live_plot_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"w");

    if( fp == NULL ){
        printf(
          "\n\n ERROR LivePlot_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{
        fclose(fp);
        return;
    }
}

void LivePlot_append(
  const char *ode_name, char *eos_name, double *pars, double R, int n
){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/live_plot_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"a");

    if( fp == NULL ){
        printf(
          "\n\n ERROR LivePlot_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{

        if(kmax){
            fprintf(
              fp,
              "# phiScal_c = %e, R = %e \n",
              yp[1][1], R
            );

            for(int i=1; i <= kmax && xp[i]; i++){
                fprintf(fp,"%e ", xp[i]);
                xp[i] = 0;
                for(int k=1; k <= n; k++){
                    fprintf(fp,"%e ", yp[k][i]);
                    yp[k][i] = 0;
                }

                fprintf(fp,"%e \n", rhop[i]);
                rhop[i] = 0;
            }

        }else{
            fprintf(fp,"# NO POINTS OOOO");
        }

        fclose(fp);
    }
}

void ResultFile_open(const char *ode_name, char *eos_name, double *pars){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/STT_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"w");

    if( fp == NULL ){
        printf(
          "\n\n ERROR ResultFile_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{

        fprintf(fp,"# p_c phiScal_c M AR rho_c delta_phiScal \n");

        fclose(fp);
        return;
    }
}

void ResultFile_append(
  const char *ode_name, char *eos_name, double *pars,
  double p_c, double phiScal_c, double M, double AR, double rho_c,
  double delta_phiScal
){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/STT_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"a");

    if( fp == NULL ){
        printf(
          "\n\n ERROR LivePlot_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{

        fprintf(
          fp,
          "%e %e %e %e %e %e\n",
          p_c, phiScal_c, M, AR, rho_c, delta_phiScal
        );

        fclose(fp);
    }
}

void LivePlot_phiScal_J_open(const char *ode_name, char *eos_name, double *pars){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/live_plot_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"w");

    if( fp == NULL ){
        printf(
          "\n\n ERROR LivePlot_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{
        fclose(fp);
        return;
    }
}

void LivePlot_phiScal_J_append(
  const char *ode_name, char *eos_name, double *pars, double R, int n
){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/live_plot_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"a");

    if( fp == NULL ){
        printf(
          "\n\n ERROR LivePlot_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{

        if(kmax){
            fprintf(
              fp,
              "# phiScal_c = %e, R = %e \n",
              yp[1][1], R
            );

            for(int i=1; i <= kmax && xp[i]; i++){
                fprintf(fp,"%e ", xp[i]);
                xp[i] = 0;
                for(int k=1; k <= n; k++){
                    fprintf(fp,"%e ", yp[k][i]);
                    yp[k][i] = 0;
                }

                fprintf(fp,"%e \n", rhop[i]);
                rhop[i] = 0;
            }

        }else{
            fprintf(fp,"# NO POINTS OOOO");
        }

        fclose(fp);
    }
}

void ResultFile_phiScal_J_open(const char *ode_name, char *eos_name, double *pars){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/STT_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"w");

    if( fp == NULL ){
        printf(
          "\n\n ERROR ResultFile_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{

        fprintf(fp,"# p_c phiScal_c M AR rho_c J \n");

        fclose(fp);
        return;
    }
}

void ResultFile_phiScal_J_append(
  const char *ode_name, char *eos_name, double *pars,
  double p_c, double phiScal_c, double M, double AR, double rho_c, double J
){

    char \
      path[256] = "/home/dimitar/projects/STT_theories/results/STT_",
      beta[16], m[16], lambda[16];

    strcat(path, ode_name);
    strcat(path, "_");
    strcat(path, eos_name);

    strcat(path, "_beta");
    sprintf( beta,"%.3e", pars[1] );
    strcat(path, beta);

    strcat(path, "_m");
    sprintf( m,"%.3e", pars[2] );
    strcat(path, m);

    strcat(path, "_lambda");
    sprintf( lambda,"%.3e", pars[3] );
    strcat(path, lambda);

    FILE *fp = fopen(path,"a");

    if( fp == NULL ){
        printf(
          "\n\n ERROR LivePlot_open oppenning file:\n\n \
          %s ERROR line 178\n\n",
          path
        );
        perror("Error");

        exit(178);
    }else{

        fprintf(
          fp,
          "%e %e %e %e %e %e \n",
          p_c, phiScal_c, M, AR, rho_c, J
        );

        fclose(fp);
    }
}

#undef PARAM_VAL_CHAR_LENG
#undef DEBUGGING_dvector_copy
#undef DEBUGGING_OPEN_WRITE_RESULT
#undef DEBUGGING_dvector_copy_to_index
