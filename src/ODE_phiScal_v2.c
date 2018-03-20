#include "ExternalHeaders.h"

#define COUNT_EQS 5
#define ODE_NAME "phiScal"

// interval of central pressures to go over through
#define P_START 1e-5 // 1e-5
#define P_END 5e-3 // 5e-3

static double \
  // values to save in the ResultFile
  AR, R, rho_c,
  // the parameters of the system
  *GV_PARAMETERS_VALUES,
  // we use this variable to set the central pressure in initial_values_init
  p_current = 3e-4,
  // tiny parameter for our working pleasure
  tiny = 1e-30;

static EOSmodelInfoStruct *eos;

static int minimal_p_power = -16;

static double get_power(double x){

    return floor(log10(fabs( x + tiny)));

}

static void phiScal_foo(double x, double *y, double *dydx){

    double \
        beta_phiScal = GV_PARAMETERS_VALUES[1], \
        m_phiScal = GV_PARAMETERS_VALUES[2], \
        lambda_phiScal = GV_PARAMETERS_VALUES[3], \
        r = x, \
        phiScal = y[1], \
        Q = y[2], \
        p = y[3], \
        LambdaMetr = y[4], \
        m = y[5], \
        \
        Vhat = 2*pow(m_phiScal,2)*pow(phiScal,2) + lambda_phiScal*pow(phiScal,4), \
        Vhat_dphiScal = 4*pow(m_phiScal,2)*phiScal + 4*lambda_phiScal*pow(phiScal,3), \
        alpha = beta_phiScal*phiScal, \
        A = exp(1.0/2 * beta_phiScal * pow(phiScal,2)), \
        exp_2LambdaMetr = exp(2*LambdaMetr ), \
        \
        step4_A = pow(A,4), \
        rho;

    if( !R && p &&  get_power(p) <= minimal_p_power ){
        R = r;
        AR = A*r;
        p = 0;
        rho = 0;
    }else if( get_power(p) <= minimal_p_power ){
        p = 0;
        rho = 0;
    }else{
        EOSeq(eos,p);
        rho = eos->current;

        if(!rho_c){
            rho_c = rho;
        }
    }

    rho_tmp = rho;

    double \
      PhiMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*p*exp_2LambdaMetr
          + pow(Q,2)
          - 1.0/2*exp_2LambdaMetr*Vhat
          + pow(1/r,2)*(exp_2LambdaMetr-1)
        ), \
      phiScal_dr = Q, \
      LambdaMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*rho*exp_2LambdaMetr
          + pow(Q,2)
          + 1.0/2*exp_2LambdaMetr*Vhat
          - pow(1/r,2)*(exp_2LambdaMetr-1)
        ), \
      Q_dr = \
        4*GV_PI*alpha*step4_A*(rho - 3*p)*exp_2LambdaMetr \
        + 1.0/4*Vhat_dphiScal*exp_2LambdaMetr \
        - (PhiMetr_dr - LambdaMetr_dr + 2/r )*Q, \
      p_dr = -( rho + p ) * ( PhiMetr_dr + alpha*Q ), \
      m_dr = \
        pow(r,2) * (
          4 * GV_PI * step4_A* rho
          + 1.0/2 * 1/exp_2LambdaMetr * pow(Q,2)
          + 1.0/4 * Vhat
        );

    dydx[1] = phiScal_dr;
    dydx[2] = Q_dr;
    dydx[3] = p_dr;
    dydx[4] = LambdaMetr_dr;
    dydx[5] = m_dr;

    return;
}

static void phiScal_parameters_init(void){

    GV_PARAMETERS_VALUES = dvector(1,3);

    // beta
    GV_PARAMETERS_VALUES[1] = -6;

    // m
    GV_PARAMETERS_VALUES[2] = 0;

    // lambda
    GV_PARAMETERS_VALUES[3] = 0;

    return;

}

static void odeint_point_arrs_init(void){

    // how many points we want to be saved, if any
    // creates the arrays which are used by odeint
    // rhop is for the density, as it is not part of the system itself
    kmax = 1e4;

    if(kmax){
        xp = dvector(1, kmax);
        yp = dmatrix(1, COUNT_EQS, 1, kmax);
        rhop = dvector(1, kmax);
        dxsav = 1e-1;
    }

    return;
}

static void odeint_point_arrs_free(void){

    // how many points we want to be saved, if any
    // creates the arrays which are used by odeint
    // rhop is for the density, as it is not part of the system itself
    if(kmax){
        free(xp);
        free_dmatrix( yp, 1, COUNT_EQS, 1, kmax );
        free(rhop);
    }

    return;
}

static double *initial_values_init(void){

    double *y0 = dvector(1,COUNT_EQS);

    // phiScal
    y0[1] = 0;

    // Q = dphiScaldr
    y0[2] = 0;

    // pressure at the center
    y0[3] = p_current;

    // Lambda - metric function
    y0[4] = 0;

    // the mass of the star
    y0[5] = 0;

    return y0;
}

static ODE_struct *ODE_struct_init(void){

    ODE_struct *_ode = calloc(1, sizeof(ODE_struct));

    // set the initial values of the system
    _ode->y = initial_values_init();

    // start of the integration interval
    _ode->x_start = 1e-30;

    // end of the integration interval
    _ode->x_end = 1e4;

    return _ode;
}

static void ODE_struct_free( ODE_struct **_ode ){

    free((*_ode)->y);
    free((*_ode));

    return;
}

void single_integrate_phiScal(void){

    ODE_struct *_ode = ODE_struct_init();

    // dummy variables to keep track for the amount of ok and bad steps
    // done by odeint
    int nok = 0, nbad = 0;

    // init the GV_PARAMETERS_VALUES for the scalar field
    // beta, m, lambda
    phiScal_parameters_init();

    // init the xp, yp, kmax, dxsave for the odeint
    // they are used for saving kmax values of the system at relative dxsave steps
    odeint_point_arrs_init();

    // init the equation of state
    eos = calloc(1,sizeof(EOSmodelInfoStruct));
    eos_init(&eos);

    // global variables, they are to keep the radii of the star
    // and the central density
    R = 0;
    AR = 0;
    rho_c = 0;

    // global variables defined in odeint
    // it will tell us if the system went boom and where
    did_it_go_boom = 0;
    where_it_went_boom = 0;

    printf("\n starting at %e \n", _ode->x_start);
    for(int i = 1; i <= COUNT_EQS; i++){
        printf("\t [ %d ] = %e \n", i, _ode->y[i]);
    }

    LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    odeint(
      _ode->y, COUNT_EQS,
      _ode->x_start, _ode->x_end,
      &nok, &nbad,
      &phiScal_foo
    );

    printf("\n end at %e \n", _ode->x_end);
    for(int i = 1; i <= COUNT_EQS; i++){
        printf("\t [ %d ] = %e \n", i, _ode->y[i]);
    }

    LivePlot_append(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES, R, COUNT_EQS);

    ODE_struct_free(&_ode);
    free(GV_PARAMETERS_VALUES);
    odeint_point_arrs_free();
    eos_free(&eos);

    return;
}

static void integrate_phiScal( ODE_struct *_ode ){

    int nok = 0, nbad = 0;

    // init the xp, yp, kmax, dxsave for the odeint
    // they are used for saving kmax values of the system at relative dxsave steps
    odeint_point_arrs_init();

    // global variables, they are to keep the radii of the star
    // and the central density
    R = 0;
    AR = 0;
    rho_c = 0;

    // global variables defined in odeint
    // it will tell us if the system went boom and where
    did_it_go_boom = 0;
    where_it_went_boom = 0;

    odeint(
      _ode->y, COUNT_EQS,
      _ode->x_start, _ode->x_end,
      &nok, &nbad,
      &phiScal_foo
    );

    // let us know if the system went boom
    if(did_it_go_boom){
        printf(
          "\n integrate_phiscal says system went boom at %e \n",
          where_it_went_boom
        );

        LivePlot_append(
          ODE_NAME,
          eos->model_name,
          GV_PARAMETERS_VALUES,
          R,
          COUNT_EQS
        );

        exit(123);
    }

    // if this is part of fadj in solver_newt do not append it to the LivePlot
    if(!GV_fadj_switch){
        LivePlot_append(
          ODE_NAME,
          eos->model_name,
          GV_PARAMETERS_VALUES,
          R,
          COUNT_EQS
        );
    }

    odeint_point_arrs_free();

    return;
}

static void shoot_regular_execute(int n, double *v, double *f){

    // set the initial data for the ode system
    ODE_struct *_ode = ODE_struct_init();

    _ode->y[1] = v[1];

    integrate_phiScal(_ode);

    // on the left side we want scalar field to be 0
    f[1] = 0 - _ode->y[1];

    ODE_struct_free(&_ode);

    return;
}

void single_shoot_regular(void){

    // init the GV_PARAMETERS_VALUES for the scalar field
    // beta, m, lambda
    phiScal_parameters_init();

    // init the equation of state
    eos = calloc(1,sizeof(EOSmodelInfoStruct));
    eos_init(&eos);

    // reset the LivePlot file
    LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    int \
      newt_n = 1, // the size of the guess values vector v
      newt_check = 0; // variable to check if we have stumbled in a minimum

    // the guess vector itself with initial guesses
    double *newt_v = dvector(1, newt_n);

    // initial guess for the scalar field
    newt_v[1] = -5e-1;

    printf("\n v init: \n");
    for(int i=1; i <= newt_n; i++){
        printf("\t v[%d] = %e \n", i, newt_v[i]);
    }

    newt(
        newt_v,
        newt_n,
        &newt_check,
        &shoot_regular_execute
    );

    printf("\n v final: \n");
    for(int i=1; i <= newt_n; i++){
        printf("\t v[%d] = %e \n", i, newt_v[i]);
    }

    free(newt_v);
    eos_free(&eos);
    free(GV_PARAMETERS_VALUES);

    return;
}

void iterate_pressure_shoot_regular(void){

    // init the GV_PARAMETERS_VALUES for the scalar field
    // beta, m, lambda
    phiScal_parameters_init();

    // init the equation of state
    eos = calloc(1,sizeof(EOSmodelInfoStruct));
    eos_init(&eos);

    // prepare the result file for feeding information
    ResultFile_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    int \
      iterate = 1, // control if we have went over all central pressures
      newt_check = 0, // variable to check if we have stumbled at a minimum
      newt_n = 1; // the size of the guess size vector

    double \
      // central pressure values
      pressure_start = P_START, pressure_end = P_END,
      // the guess vector itself
      *newt_v = dvector(1,newt_n),
      // dummy phiScal var
      phiScal_c;

    p_current = pressure_start;
    newt_v[1] = -5e-2;
    phiScal_c = newt_v[1];

    while(iterate){

        newt_v[1] = phiScal_c;

        printf("\n p_c = %.3e \n", p_current );
        for(int i=1; i <= newt_n; i++){
            printf("\t v[%d] = %e \n", i, newt_v[i]);
        }

        // reset the LivePlot file
        LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

        newt(
            newt_v,
            newt_n,
            &newt_check,
            &shoot_regular_execute
        );

        printf("\n\t\t one final integration to check boundary conditions \n");

        ODE_struct *_ode = ODE_struct_init();
        _ode->y[1] = newt_v[1];

        integrate_phiScal(_ode);

        printf("\n v final (newt check -> %d): \n", newt_check);

        printf(
          "\t v[%d] = %e --> delta = %e \n", 0, newt_v[1], 0 - _ode->y[1]
        );

        ResultFile_append(
          ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES,
          p_current, newt_v[1], _ode->y[5], AR, rho_c,
          0 - _ode->y[1]
        );

        ODE_struct_free(&_ode);

        if(fabs(newt_v[1]) > fabs(phiScal_c)){
            printf(
              "\n\t\t changing the initial guess from %e to %e \n",
              phiScal_c,
              newt_v[1]
            );

            phiScal_c = newt_v[1];

        }else if( get_power(newt_v[1] + tiny) >= -2 ){
            printf(
              "\n\t\t changing the initial guess from %e to %e \n",
              phiScal_c,
              newt_v[1]
            );

            phiScal_c = newt_v[1];
        }

        int current_power = (int)get_power(p_current);
        switch(current_power){
            case -5:
                p_current += 2*pow10( current_power - 1 );
                break;
            case -4:
                p_current += 1*pow10( current_power - 1 );
                break;
            case -3:
                p_current += 5*pow10( current_power - 2 );
                break;
            default:
                printf("\n %e not known, terminating... \n", p_current);
                exit(1376);
                break;
        }

        if(p_current > pressure_end){
            iterate = 0;
        }else{
            iterate = 1;
        }

    }

    free(newt_v);
    eos_free(&eos);
    free(GV_PARAMETERS_VALUES);

    return;
}

