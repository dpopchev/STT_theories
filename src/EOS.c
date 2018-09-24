#include "ExternalHeaders.h"

#define DEBUGGING_EOS_INIT 0

// the total amount of eos we have
#define DEF_TOTAL_EOS_NAME_COUNT 35

// which model we will use
#define DEF_EOS_MODEL_NUM 3
// if we will use splines, the dimensions of the vector
#define DEF_EOS_PICEWISE_VECTORS_DIM 7

const char \
    *default_eos_names[] = {
        "", "EOSII",
        "PAL6",     "SLy",
        "APR1",     "APR2",      "APR3",     "APR4",
        "FPS",     "WFF1",     "WFF2",      "WFF3",
        "BBB2",     "BPAL12",     "ENG",     "MPA1",
        "MS1",      "MS2",     "MS1b",      "PS",
        "GS1",     "GS2",       "BGN1H1",     "GNH3",
        "H1",       "H2",       "H3",       "H4",
        "H5",       "H6",       "H7",       "PCL2",
        "ALF1",       "ALF2",     "ALF3",     "ALF4"
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

    //~ if(
      //~ (*eos)->model_num == 2 || (*eos)->model_num == 4  || (*eos)->model_num > 4
    //~ ){
        //~ printf(
          //~ "\n eos number %d not working yet, eos.c eos_init, exitting line 127 \n",
          //~ (*eos)->model_num
        //~ );

        //~ exit(127);
    //~ }

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

        case 2:
            tmp_log_p1 = 34.380;
            tmp_Gamma1 = 2.227;
            tmp_Gamma2 = 2.189;
            tmp_Gamma3 = 2.159;

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

        case 3:
            tmp_log_p1 = 34.384;
            tmp_Gamma1 = 3.005;
            tmp_Gamma2 = 2.988;
            tmp_Gamma3 = 2.851;

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

        case 4:
            tmp_log_p1 = 33.943;
            tmp_Gamma1 = 2.442;
            tmp_Gamma2 = 3.256;
            tmp_Gamma3 = 2.908;

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

        case 5:
            tmp_log_p1 = 34.126;
            tmp_Gamma1 = 2.643;
            tmp_Gamma2 = 3.014;
            tmp_Gamma3 = 2.945;

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

        case 6:
            tmp_log_p1 = 34.392;
            tmp_Gamma1 = 3.166;
            tmp_Gamma2 = 3.573;
            tmp_Gamma3 = 3.281;

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

        case 7:
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

        case 8:
            tmp_log_p1 = 34.283;
            tmp_Gamma1 = 2.985;
            tmp_Gamma2 = 2.863;
            tmp_Gamma3 = 2.600;

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

        case 9:
            tmp_log_p1 = 34.031;
            tmp_Gamma1 = 2.519;
            tmp_Gamma2 = 3.791;
            tmp_Gamma3 = 3.660;

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

        case 10:
            tmp_log_p1 = 34.233;
            tmp_Gamma1 = 2.888;
            tmp_Gamma2 = 3.475;
            tmp_Gamma3 = 3.517;

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

        case 11:
            tmp_log_p1 = 34.283;
            tmp_Gamma1 = 3.329;
            tmp_Gamma2 = 2.952;
            tmp_Gamma3 = 2.589;

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

        case 12:
            tmp_log_p1 = 34.331;
            tmp_Gamma1 = 3.418;
            tmp_Gamma2 = 2.835;
            tmp_Gamma3 = 2.832;

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

        case 13:
            tmp_log_p1 = 34.358;
            tmp_Gamma1 = 2.209;
            tmp_Gamma2 = 2.201;
            tmp_Gamma3 = 2.176;

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

        case 14:
            tmp_log_p1 = 34.437;
            tmp_Gamma1 = 3.514;
            tmp_Gamma2 = 3.130;
            tmp_Gamma3 = 3.168;

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

        case 15:
            tmp_log_p1 = 34.495;
            tmp_Gamma1 = 3.446;
            tmp_Gamma2 = 3.572;
            tmp_Gamma3 = 2.887;

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

        case 16:
            tmp_log_p1 = 34.858;
            tmp_Gamma1 = 3.224;
            tmp_Gamma2 = 3.033;
            tmp_Gamma3 = 1.325;

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

        case 17:
            tmp_log_p1 = 34.605;
            tmp_Gamma1 = 2.447;
            tmp_Gamma2 = 2.184;
            tmp_Gamma3 = 1.855;

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

        case 18:
            tmp_log_p1 = 34.855;
            tmp_Gamma1 = 3.456;
            tmp_Gamma2 = 3.011;
            tmp_Gamma3 = 1.425;

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

        case 19:
            tmp_log_p1 = 34.671;
            tmp_Gamma1 = 2.216;
            tmp_Gamma2 = 1.640;
            tmp_Gamma3 = 2.365;

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

        case 20:
            tmp_log_p1 = 34.504;
            tmp_Gamma1 = 2.350;
            tmp_Gamma2 = 1.267;
            tmp_Gamma3 = 2.421;

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

        case 21:
            tmp_log_p1 = 34.642;
            tmp_Gamma1 = 2.519;
            tmp_Gamma2 = 1.571;
            tmp_Gamma3 = 2.314;

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

        case 22:
            tmp_log_p1 = 34.623;
            tmp_Gamma1 = 3.258;
            tmp_Gamma2 = 1.472;
            tmp_Gamma3 = 2.464;

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

        case 23:
            tmp_log_p1 = 34.648;
            tmp_Gamma1 = 2.664;
            tmp_Gamma2 = 2.194;
            tmp_Gamma3 = 2.304;

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

        case 24:
            tmp_log_p1 = 34.564;
            tmp_Gamma1 = 2.595;
            tmp_Gamma2 = 1.845;
            tmp_Gamma3 = 1.897;

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

        case 25:
            tmp_log_p1 = 34.617;
            tmp_Gamma1 = 2.775;
            tmp_Gamma2 = 1.855;
            tmp_Gamma3 = 1.858;

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

        case 26:
            tmp_log_p1 = 34.646;
            tmp_Gamma1 = 2.787;
            tmp_Gamma2 = 1.951;
            tmp_Gamma3 = 1.901;

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

        case 27:
            tmp_log_p1 = 34.669;
            tmp_Gamma1 = 2.909;
            tmp_Gamma2 = 2.246;
            tmp_Gamma3 = 2.144;

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

        case 28:
            tmp_log_p1 = 34.609;
            tmp_Gamma1 = 2.793;
            tmp_Gamma2 = 1.974;
            tmp_Gamma3 = 1.915;

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

        case 29:
            tmp_log_p1 = 34.593;
            tmp_Gamma1 = 2.637;
            tmp_Gamma2 = 2.121;
            tmp_Gamma3 = 2.064;

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

        case 30:
            tmp_log_p1 = 34.559;
            tmp_Gamma1 = 2.621;
            tmp_Gamma2 = 2.048;
            tmp_Gamma3 = 2.006;

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

        case 31:
            tmp_log_p1 = 34.507;
            tmp_Gamma1 = 2.554;
            tmp_Gamma2 = 1.880;
            tmp_Gamma3 = 1.977;

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

        case 32:
            tmp_log_p1 = 34.055;
            tmp_Gamma1 = 2.013;
            tmp_Gamma2 = 3.389;
            tmp_Gamma3 = 2.033;

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

        case 33:
            tmp_log_p1 = 34.616;
            tmp_Gamma1 = 4.070;
            tmp_Gamma2 = 2.411;
            tmp_Gamma3 = 1.890;

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

        case 34:
            tmp_log_p1 = 34.283;
            tmp_Gamma1 = 2.883;
            tmp_Gamma2 = 2.653;
            tmp_Gamma3 = 1.952;

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

        case 35:
            tmp_log_p1 = 34.314;
            tmp_Gamma1 = 3.009;
            tmp_Gamma2 = 3.438;
            tmp_Gamma3 = 1.803;

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

    if(eos->model_num >= 2){
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
    }else{
        printf("\n EOSeq do not know what to do %d \n", eos->model_num);

        return;
    }

    //~ switch(eos->model_num){

        //~ case 1:{

            //~ double p0 = 2.7039e-4, gamma = 2.34, k = 0.0195;

            //~ eos->current = p0*pow( pressure/(k*p0), 1.0/gamma);

            //~ return;
            //~ break;

        //~ }default:
            //~ printf(
              //~ "\n EOS.c eos_init %d unknown, exiting... line 260\n", eos->model_num
            //~ );


            //~ exit(260);
            //~ break;
    //~ }

    return;
}

#undef DEBUGGING_EOS_INIT
#undef DEF_TOTAL_EOS_NAME_COUNT
#undef DEF_EOS_MODEL_NUM
#undef DEF_EOS_PICEWISE_VECTORS_DIM
