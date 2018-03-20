#include "ExternalHeaders.h"

#define DEBUGGING_EOS_INIT 0

// the total amount of eos we have
#define DEF_TOTAL_EOS_NAME_COUNT 4

// which model we will use
#define DEF_EOS_MODEL_NUM 3

// if we will use splines, the dimensions of the vector
#define DEF_EOS_PICEWISE_VECTORS_DIM 7

const char \
    *default_eos_names[] = {
      "", "EOSII", "Sly4", "AkmalPR", "piecewAPR4"
    };

static void picewise_EOS_load(
  double log_p1, double Gamma1, double Gamma2,
  double Gamma3, EOSmodelInfoStruct *EOSdata
){

    double *K, *Gamma, *a, *divid_rho, *divid_p, *N, LV_C = 2.9979e10;

    K = dvector(0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    Gamma = dvector(0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    a = dvector(0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    divid_rho = dvector(0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    divid_p = dvector(0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    N = dvector(0,DEF_EOS_PICEWISE_VECTORS_DIM-1);

    a[0] = 0;

    K[0] = 6.80110e-9;
    Gamma[0] = 1.58425;
    divid_rho[0] = 2.44034e7;

    K[1] = 1.06186e-6;
    Gamma[1] = 1.28733;
    divid_rho[1] = 3.78358e11;

    K[2] = 5.32697e+1;
    Gamma[2] = 0.62223;
    divid_rho[2] = 2.62780e12;

    K[3] = 3.99874e-8;
    Gamma[3] = 1.35692;
    divid_rho[3] = 0.00000e00;

    for(int i = 0; i <= 3; i++){
        K[i] *= pow(LV_C, 2.0);
    }

    K[4] = 0.00000e00;
    Gamma[4] = Gamma1;
    divid_rho[4] = pow10(14.7);

    K[5] = 0.00000e00;
    Gamma[5] = Gamma2;
    divid_rho[5] = 1.00000e15;

    K[6] = 0.00000e00;
    Gamma[6] = Gamma3;
    divid_rho[6] = 9.99999e99;

    K[4] = pow10(log_p1) / pow(divid_rho[4], Gamma[4]);

    if( Gamma[4] == Gamma[3] ){
        divid_rho[3] = 1e14;
    }
    else{
        divid_rho[3] = pow( K[3] / K[4], 1.0 / (Gamma[4] - Gamma[3]) );
    }

    for(int i = 5; i <= 6; i++){
        K[i] = K[i-1]*pow(divid_rho[i-1], Gamma[i-1] - Gamma[i]);
    }

    for(int i = 1; i <= 6; i++){
        a[i] = \
          a[i-1] \
          + K[i-1]/pow(LV_C, 2.0)/(Gamma[i-1] - 1) \
          *pow(divid_rho[i-1], Gamma[i-1] - 1) \
          - K[i]/pow(LV_C, 2.0)/(Gamma[i] - 1) \
          *pow(divid_rho[i-1], Gamma[i]- 1);
    }


    for(int i = 0; i <= 6; i++){
        divid_p[i] = K[i]*pow(divid_rho[i],Gamma[i]);
    }

    for(int i=0; i<= 6; i++){
        EOSdata->PieceWise_K[i+1] = K[i];
        EOSdata->PieceWise_Gamma[i+1] = Gamma[i];
        EOSdata->PieceWise_a[i+1] = a[i];
        EOSdata->PieceWise_density[i+1] = divid_rho[i];
        EOSdata->PieceWise_presure[i+1] = divid_p[i];
        EOSdata->PieceWise_N[i+1] = N[i];
    }

    free_dvector(K,0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    free_dvector(Gamma,0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    free_dvector(a,0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    free_dvector(divid_rho,0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    free_dvector(divid_p,0,DEF_EOS_PICEWISE_VECTORS_DIM-1);
    free_dvector(N,0,DEF_EOS_PICEWISE_VECTORS_DIM-1);

    return;
}

void eos_init(EOSmodelInfoStruct **eos){

    (*eos)->model_num = DEF_EOS_MODEL_NUM;

    if(
      (*eos)->model_num == 2 || (*eos)->model_num == 4  || (*eos)->model_num > 4
    ){
        printf(
          "\n eos number %d not working yet, eos.c eos_init, exitting line 127 \n",
          (*eos)->model_num
        );

        exit(127);
    }

    (*eos)->model_name = calloc(32,sizeof(char));
    strcpy( (*eos)->model_name, default_eos_names[(*eos)->model_num]);

    double tmp_log_p1, tmp_Gamma1, tmp_Gamma2, tmp_Gamma3;

    switch((*eos)->model_num){

        case 1:

            (*eos)->PieceWise_K = NULL;
            (*eos)->PieceWise_Gamma = NULL;
            (*eos)->PieceWise_a = NULL;
            (*eos)->PieceWise_density = NULL;
            (*eos)->PieceWise_presure = NULL;
            (*eos)->PieceWise_N = NULL;

            break;

        case 3:

            tmp_log_p1 = 34.269;
            tmp_Gamma1 = 2.830;
            tmp_Gamma2 = 3.445;
            tmp_Gamma3 = 3.348;

            (*eos)->PieceWise_K = dvector(1,DEF_EOS_PICEWISE_VECTORS_DIM);
            (*eos)->PieceWise_Gamma = dvector(1,DEF_EOS_PICEWISE_VECTORS_DIM);
            (*eos)->PieceWise_a = dvector(1,DEF_EOS_PICEWISE_VECTORS_DIM);
            (*eos)->PieceWise_density = dvector(1,DEF_EOS_PICEWISE_VECTORS_DIM);
            (*eos)->PieceWise_presure = dvector(1,DEF_EOS_PICEWISE_VECTORS_DIM);
            (*eos)->PieceWise_N = dvector(1,DEF_EOS_PICEWISE_VECTORS_DIM);

            picewise_EOS_load(
              tmp_log_p1, tmp_Gamma1, tmp_Gamma2, tmp_Gamma3, (*eos)
            );

            break;

        default:
            printf(
              "\n EOS.c eos_init %d unknown, exiting... line 171\n", (*eos)->model_num
            );


            exit(171);
            break;
    }

    return;
}

void eos_print_info_stdout(EOSmodelInfoStruct *eos){

    printf(
        "\n EOS info: \n"
        "\n\t General info: \n"
        "\t\t name %s \n",
        eos->model_name
    );

    return;

}

void eos_print_info_ResultFile( EOSmodelInfoStruct *eos, FILE *fp ){

    fprintf(
        fp,
        "\n EOS info: \n"
        "\n\t General info: \n"
        "\t\t name %s \n",
        eos->model_name
    );

    fclose(fp);

    return;
}

void eos_free(EOSmodelInfoStruct **eos){

    free((*eos)->model_name);
    free((*eos)->PieceWise_K);
    free((*eos)->PieceWise_Gamma);
    free((*eos)->PieceWise_a);
    free((*eos)->PieceWise_density);
    free((*eos)->PieceWise_presure);
    free((*eos)->PieceWise_N);
    free(*eos);

    return;
}

void EOSeq( EOSmodelInfoStruct *eos, double pressure ){

    switch(eos->model_num){

        case 1:{

            double p0 = 2.7039e-4, gamma = 2.34, k = 0.0195;

            eos->current = p0*pow( pressure/(k*p0), 1.0/gamma);

            return;
            break;

        }case 3:

            pressure *= GV_PRESSURE_UNITS;
            int i = 1;
            double LV_C = 2.9979e10;

            while( pressure > eos->PieceWise_presure[i] ){
                i++;
            }

            double eps =  \
              ( 1.0 + eos->PieceWise_a[i] ) \
              * pow( pressure / eos->PieceWise_K[i],1.0 / eos->PieceWise_Gamma[i])\
              + pressure/pow( LV_C, 2 ) / ( eos->PieceWise_Gamma[i] - 1.0 );

            eos->current = eps/GV_DENSITY_UNITS;

            return;
            break;

        default:
            printf(
              "\n EOS.c eos_init %d unknown, exiting... line 260\n", eos->model_num
            );


            exit(260);
            break;
    }

    return;
}

#undef DEBUGGING_EOS_INIT
#undef DEF_TOTAL_EOS_NAME_COUNT
#undef DEF_EOS_MODEL_NUM
#undef DEF_EOS_PICEWISE_VECTORS_DIM
