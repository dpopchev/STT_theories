#include "ExternalHeaders.h"

#define COUNT_EQS 5
#define ODE_NAME "phiScal"

// interval of central pressures to go over through
#define P_START 2.2e-4 // 1e-5
#define P_END 2e-3 // 5e-3

// value of the infinity to use
#define R_INF 50

// parameters of the scalar field
#define BETA -6
#define M 2e-2
#define LAMBDA 0

static double \
  // values to save in the ResultFile
  AR, R, rho_c,
  // the parameters of the system
  *GV_PARAMETERS_VALUES,
  // we use this variable to set the central pressure in initial_values_init
  p_current = P_START,
  // set where the infinity should be
  r_inf = R_INF,
  // tiny parameter for our working pleasure
  tiny = 1e-30,
  // switch to control if we are just getting the central scalar field
  // if we do, we dont want to save anything
  GV_GET_PHISCAL = 0;

static EOSmodelInfoStruct *eos;

static int minimal_p_power = -20;

// returns the power of the number n
static double get_power(double x){

    return floor(log10(fabs( x + tiny)));

}

// system to integrate - phiScal + mass
static void phiScal_foo(double x, double *y, double *dydx){

    double \
        beta_phiScal = GV_PARAMETERS_VALUES[1],
        m_phiScal = GV_PARAMETERS_VALUES[2],
        lambda_phiScal = GV_PARAMETERS_VALUES[3],
        r = x,
        phiScal = y[1],
        Q = y[2],
        p = y[3],
        LambdaMetr = y[4],
        m = y[5],
        \
        Vhat = 2*pow(m_phiScal,2)*pow(phiScal,2) + lambda_phiScal*pow(phiScal,4),
        Vhat_dphiScal = 4*pow(m_phiScal,2)*phiScal + 4*lambda_phiScal*pow(phiScal,3),
        alpha = beta_phiScal*phiScal,
        A = exp(1.0/2 * beta_phiScal * pow(phiScal,2)),
        exp_2LambdaMetr = exp(2*LambdaMetr ),
        \
        step4_A = pow(A,4),
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

// set the parameters of the scalar field here
static void phiScal_parameters_init(void){

    GV_PARAMETERS_VALUES = dvector(1,3);

    // beta
    GV_PARAMETERS_VALUES[1] = BETA;

    // m
    GV_PARAMETERS_VALUES[2] = M;

    // lambda
    GV_PARAMETERS_VALUES[3] = LAMBDA;

    return;

}

// initializes the odeint arrays xp, yp, rhop in which it saves
// the values at approximately dxsave steps
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

// frees the arrays
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

// set the initial values of the system
// used by all functions below
// when used for shooting we just substitute by hand
// the central pressure is set by the static global var p_current
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

// initializes the variables needed for odeint
// sets the start and end of the interval of integration interest
// the x_end is set by the static global variable r_inf
static ODE_struct *ODE_struct_init(void){

    ODE_struct *_ode = calloc(1, sizeof(ODE_struct));

    // set the initial values of the system
    _ode->y = initial_values_init();

    // start of the integration interval
    _ode->x_start = 1e-30;

    // end of the integration interval
    _ode->x_end = r_inf;

    return _ode;
}

// free the ODE_struct
static void ODE_struct_free( ODE_struct **_ode ){

    free((*_ode)->y);
    free((*_ode));

    return;
}

// integrate just ones with provided initial values
// and integration interval
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

// integration function used by the shooting functions below
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

// shoot one time with provided initial guess
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
    newt_v[1] = -5e-2;

    printf("\n p_c = %.3e v init: \n", p_current);
    for(int i=1; i <= newt_n; i++){
        printf("\t v[%d] = %e \n", i, newt_v[i]);
    }

    newt(
        newt_v,
        newt_n,
        &newt_check,
        &shoot_regular_execute
    );

    printf("\n\t v final: \n");
    for(int i=1; i <= newt_n; i++){
        printf("\t v[%d] = %e \n", i, newt_v[i]);
    }

    free(newt_v);
    eos_free(&eos);
    free(GV_PARAMETERS_VALUES);

    return;
}

// check if the achieved central phiScala is greater then the previous one
// or it is of the power -2
// to know if we are gonna use it for further shooting
static void check_achieved_phiScal( double *v_phiScal, double *tmp_phiScal){

    if(fabs(*v_phiScal) > fabs(*tmp_phiScal)){
        printf(
          "\n\t\t changing the initial guess from %e to %e \n",
          *tmp_phiScal,
          *v_phiScal
        );

        *tmp_phiScal = *v_phiScal;

    }else if( get_power(*v_phiScal + tiny) >= -2 ){
        printf(
          "\n\t\t changing the initial guess from %e to %e \n",
          *tmp_phiScal,
          *v_phiScal
        );

        *tmp_phiScal = *v_phiScal;
    }

    return;
}

// here we are setting the desired central pressure power-wise
// depending of the power it is we increment it with different sized step
static void change_current_pressure(void){

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

    return;
}

// iterate over the central pressure values of interest
// defined by P_START and P_END
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
    newt_v[1] = -3.065399e-02;
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

        // stop for further investigation if the boundary condition for the
        // scalar field is not satisfied
        // should be connected to the solver_newt, but I prefer to set it here
        // as the power of -10 by hand
        if(get_power(0 - _ode->y[1]) > -9){
            printf("\n scalar field boundary condition not met, terminating...");
            exit(123);
        }

        printf("\n\t v final (newt check -> %d): \n", newt_check);

        printf(
          "\t v[%d] = %e --> delta = %e \n", 1, newt_v[1], 0 - _ode->y[1]
        );

        ResultFile_append(
          ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES,
          p_current, newt_v[1], _ode->y[5], AR, rho_c,
          0 - _ode->y[1]
        );

        // stop for further investigation if the boundary condition for the
        // scalar field is not satisfied
        // should be connected to the solver_newt, but I prefer to set it here
        // as the power of -10 by hand
        if(get_power(0 - _ode->y[1]) > -9){
            printf("\n scalar field boundary condition not met, terminating...");
            exit(123);
        }

        ODE_struct_free(&_ode);

        // if achieved central scalar field is greater then the previous one - change them
        // else if it is the first to be of the power of -2 - change it
        check_achieved_phiScal( &newt_v[1], &phiScal_c);

        // increment the central pressure power wise
        change_current_pressure();

        // if the current central pressure is greater the greatest desired one
        // just end the loop
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

// shoot one time with provided initial guess
// also increase the infinity until significant changes disappear
// use the centra phiScal and integrate to first inf
// after that set phiScal and its derivative as 0 and continue integrating
void single_shoot_regular_iterate_inf(void){

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
    double \
      *newt_v = dvector(1, newt_n),
      phiScal_c = -2.605183e-02;

    // initial guess for the scalar field
    newt_v[1] = phiScal_c;

    printf("\n p_c = %.3e \n", p_current);

    int count = 100;
    while(count--){
        printf("\n v_%d for inf %e init: \n", count, r_inf);
        for(int i=1; i <= newt_n; i++){
            printf("\t v[%d] = %e \n", i, newt_v[i]);
        }

        newt(
            newt_v,
            newt_n,
            &newt_check,
            &shoot_regular_execute
        );

        printf("\n\t v_%d final: \n", count);
        for(int i=1; i <= newt_n; i++){
            printf("\t v[%d] = %e \n", i, newt_v[i]);
        }

        if(get_power(phiScal_c - newt_v[1]) >= -6){
            r_inf += 0.1*r_inf;
            phiScal_c = newt_v[1];
        }else{
            printf("No significant difference");
            count = 0;
        }
    }

    printf("\n\t\t one final integration to check boundary conditions \n");

    LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    ODE_struct *_ode = ODE_struct_init();
    _ode->y[1] = newt_v[1];

    integrate_phiScal(_ode);

    // stop for further investigation if the boundary condition for the
    // scalar field is not satisfied
    // should be connected to the solver_newt, but I prefer to set it here
    // as the power of -10 by hand
    if(get_power(0 - _ode->y[1]) > -9){
        printf("\n scalar field boundary condition not met, terminating...");
        exit(123);
    }

    printf("\n\t v final (newt check -> %d): \n", newt_check);

    printf(
      "\t v[%d] = %e --> delta = %e \n", 1, newt_v[1], 0 - _ode->y[1]
    );

    _ode->y[1] = 0;
    _ode->y[2] = 0;
    _ode->x_start = _ode->x_end;
    _ode->x_end = 1e5;

    integrate_phiScal(_ode);

    ODE_struct_free(&_ode);
    free(newt_v);
    eos_free(&eos);
    free(GV_PARAMETERS_VALUES);

    return;
}

// iterate over the central pressure values of interest
// defined by P_START and P_END
// also iterate over the r_inf until no significant difference
// occurse in the phiScal_c
void iterate_pressure_inf_shoot_regular(void){

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
    newt_v[1] = -5e-02;
    phiScal_c = newt_v[1];

    while(iterate){

        newt_v[1] = phiScal_c;

        printf("\n p_c = %.3e \n", p_current );

        // save how many incrementation of the infinity we have done
        // also control parameter for increase
        int count = 1;

        // but of course lets reset the rad coordinate infinity
        r_inf = R_INF;

        while(++count){

            double phiScal_c_tmp = newt_v[1];

            for(int i=1; i <= newt_n; i++){
                printf(
                  "\t for %e v_%d [ %d ] = %e \n",
                  r_inf, count, i, newt_v[i]
                );
            }

            // reset the LivePlot file
            LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

            newt(
              newt_v,
              newt_n,
              &newt_check,
              &shoot_regular_execute
            );

            printf("\n\t v_%d final: \n", count);
            for(int i=1; i <= newt_n; i++){
                printf("\t v[%d] = %e \n", i, newt_v[i]);
            }

            if(get_power(phiScal_c_tmp - newt_v[1]) >= -6){
                r_inf += 0.1*r_inf;
                phiScal_c_tmp = newt_v[1];
            }else{
                printf("No significant difference \n");
                break;
            }
        }

        printf("\n\t\t one final integration to check boundary conditions \n");

        ODE_struct *_ode = ODE_struct_init();
        _ode->y[1] = newt_v[1];

        integrate_phiScal(_ode);

        // stop for further investigation if the boundary condition for the
        // scalar field is not satisfied
        // should be connected to the solver_newt, but I prefer to set it here
        // as the power of -10 by hand
        if(get_power(0 - _ode->y[1]) > -9){
            printf("\n scalar field boundary condition not met, terminating...");
            exit(123);
        }

        printf("\n\t v final (newt check -> %d): \n", newt_check);

        printf(
          "\t v[%d] = %e --> delta = %e \n", 1, newt_v[1], 0 - _ode->y[1]
        );

        ResultFile_append(
          ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES,
          p_current, newt_v[1], _ode->y[5], AR, rho_c,
          0 - _ode->y[1]
        );

        // stop for further investigation if the boundary condition for the
        // scalar field is not satisfied
        // should be connected to the solver_newt, but I prefer to set it here
        // as the power of -10 by hand
        if(get_power(0 - _ode->y[1]) > -9){
            printf("\n scalar field boundary condition not met, terminating...");
            exit(123);
        }

        ODE_struct_free(&_ode);

        // if achieved central scalar field is greater then the previous one - change them
        // else if it is the first to be of the power of -2 - change it
        check_achieved_phiScal( &newt_v[1], &phiScal_c);

        // increment the central pressure power wise
        change_current_pressure();

        // if the current central pressure is greater the greatest desired one
        // just end the loop
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

// fitting point execution procedure
// if needed will stay here
static void shoot_fitting_execute(int n, double *v, double *f){

    double r_fit;

    // set the initial data left to right through fitting
    ODE_struct *_ode1 = ODE_struct_init();

    _ode1->y[1] = v[1];
    _ode1->x_end = r_fit;

    integrate_phiScal(_ode1);

    // set the initial data right to left through fitting
    ODE_struct *_ode2 = ODE_struct_init();

    _ode2->y[1] = 0;
    _ode2->y[2] = 0;
    _ode2->y[3] = 0;
    _ode2->y[4] = v[2];
    _ode2->y[5] = v[3];
    _ode2->x_end = r_fit;
    _ode2->x_start = r_inf;

    integrate_phiScal(_ode2);

    // on the left side we want scalar field to be 0
    f[1] = _ode1->y[1] - _ode2->y[1];
    f[2] = _ode1->y[4] - _ode2->y[4];
    f[3] = _ode1->y[5] - _ode2->y[5];

    ODE_struct_free(&_ode1);
    ODE_struct_free(&_ode2);

    return;
}

// the equations for slow rotation depend on those for phiScal but not vice verse
// for provided GV_PARAMETERS_VALUES, pressure and Equation of State
// override the given phiScal and inf
void get_phiScal_cVal_infVal(
  double *parm_vals, double pressure, EOSmodelInfoStruct *_eos,
  double *phiScal, double *inf, int _minimal_p_power
){

    GV_PARAMETERS_VALUES = parm_vals;

    eos = _eos;

    p_current = pressure;

    int \
      newt_n = 1, // the size of the guess values vector v
      newt_check = 0; // variable to check if we have stumbled in a minimum

    // the guess vector itself with initial guesses
    double \
      *newt_v = dvector(1, newt_n),
      phiScal_c = *phiScal,
      tmpinf = r_inf;


    // the wanted minimal power of the pressure may be different in the
    // mother function, so adopt it
    int tmp_minimal_p_power = minimal_p_power;

    minimal_p_power = _minimal_p_power;

    // initial guess for the scalar field
    newt_v[1] = phiScal_c;
    r_inf = *inf;

    // increase the infinity several times until
    // the shooted central value of the scalar field becomes
    // less and less different
    // also print how many attempts we are doing and so on
    int count = 0;
    while(++count){

        // reset the file with live plotting values
        LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

        newt(
            newt_v,
            newt_n,
            &newt_check,
            &shoot_regular_execute
        );

        printf(
          "\n\t\t v_%d = %e ----> v_%d = %e ( R = %e )\n"
          "\n\t\t count = %d with inf = %e\n\n",
          count - 1, phiScal_c, count, newt_v[1], R,
          count, r_inf
        );

        // criteria for the difference in the previous and now
        // central value of the scalar field
        if(get_power(phiScal_c - newt_v[1]) > -10 ){
            tmpinf = r_inf;
            phiScal_c = newt_v[1];
            r_inf += 0.1*r_inf;
        }else{
            printf(
                "\n\t\t No significant difference assume \n\t\t\t %e and inf %e \n",
                newt_v[1], r_inf
            );
            break;
        }
    }

    *phiScal = newt_v[1];
    *inf = r_inf;

    free(newt_v);

    // reset the minimal power of this function
    minimal_p_power = tmp_minimal_p_power;

    return;

    // I wanted to check the solution by integrating with it
    // but this proofed bad practices
    // as due the stiffens of the ODE system
    // it started to diverge
    // thus it is all commented out now, after some commit will be
    // deleted
    //ODE_struct *_ode = ODE_struct_init();
    //_ode->y[1] = newt_v[1];

    //printf(
        //"\n\t\t final integr with init %e and inf %e just to save the profiles \n",
        //_ode->y[1], _ode->x_end
    //);

    //LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    //integrate_phiScal(_ode);
    //double delta_tmp = 0 - _ode->y[1];

    //ODE_struct_free(&_ode);

    //int one_time_exec = 1;
    //while(get_power(delta_tmp) > -10 && one_time_exec){
    //while(get_power(delta_tmp) > -10){

        //r_inf = 2e2;
        //r_inf -= 0.05*r_inf;

        //newt(
            //newt_v,
            //newt_n,
            //&newt_check,
            //&shoot_regular_execute
        //);

        //printf(
            //"\n\n\t\t Reduce inf to %e to fulfil the boundary condition \n"
            //"\n\t\t lets check for %e ( R = %e )\n",
            //r_inf, newt_v[1], R
        //);

        //ODE_struct *_ode2 = ODE_struct_init();
        //_ode2->y[1] = newt_v[1];

        //LivePlot_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);
        //integrate_phiScal(_ode2);
        //delta_tmp = 0 - _ode2->y[1];
        //ODE_struct_free(&_ode2);
        ////one_time_exec = 0;
    //}

    // stop for further investigation if the boundary condition for the
    // scalar field is not satisfied
    // should be connected to the solver_newt, but I prefer to set it here
    // as the power of -10 by hand
    //if(get_power(0 - _ode->y[1]) > -10){

        //printf(
          //"\n get_phiScal_cVal_infVal scalar field boundary condition not met, terminating... \n"
        //);
        //exit(123);
    //}

    //*phiScal = newt_v[1];
    //*inf = r_inf;

    //ODE_struct_free(&_ode2);
    //free(newt_v);

    //minimal_p_power = tmp_minimal_p_power;

    //return;
}
