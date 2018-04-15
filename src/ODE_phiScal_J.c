#include "ExternalHeaders.h"

#define COUNT_EQS 9
#define ODE_NAME "phiScal_J"

// interval of central pressures to go over through
#define P_START 1e-4 // 1e-5
#define P_END 2e-3 // 5e-3

// value of the infinity to use
#define R_INF_PHISCAL 3.5e2
#define R_INF 1e9

// parameters of the scalar field
#define BETA -6
#define M 1e-3
#define LAMBDA 1

static double \
  // values to save in the ResultFile
  AR, R, rho_c,
  // the parameters of the system
  *GV_PARAMETERS_VALUES,
  // we use this variable to set the central pressure in initial_values_init
  p_current = P_START,
  // set where the infinity for phiScal and the PhiMetr
  r_inf = R_INF,
  r_inf_phiscal = R_INF_PHISCAL,
  // to use this var to set after we know the central val for phiScal
  phiScal_gv,
  // tiny parameter for our working pleasure
  tiny = 1e-30,
  GV_J, GV_M;

static EOSmodelInfoStruct *eos;

static int minimal_p_power = -17;

// returns the power of the number n
static double get_power(double x){

    return floor(log10(fabs( x + tiny)));

}

// system to integrate - phiScal + mass
static void phiScal_J_foo(double x, double *y, double *dydx){

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
        PhiMetr = y[6],
        Z = y[7],
        Omega = y[8],
        J = y[9],
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
      Omega_dr = Z,
      \
      PhiMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*p*exp_2LambdaMetr
          + pow(Q,2)
          - 1.0/2*exp_2LambdaMetr*Vhat
          + pow(1/r,2)*(exp_2LambdaMetr-1)
        ),
      \
      phiScal_dr = Q,
      \
      LambdaMetr_dr = \
        r*1.0/2 * (
          8*GV_PI*step4_A*rho*exp_2LambdaMetr
          + pow(Q,2)
          + 1.0/2*exp_2LambdaMetr*Vhat
          - pow(1/r,2)*(exp_2LambdaMetr-1)
        ),
      \
      Q_dr = \
        4*GV_PI*alpha*step4_A*(rho - 3*p)*exp_2LambdaMetr \
        + 1.0/4*Vhat_dphiScal*exp_2LambdaMetr \
        - (PhiMetr_dr - LambdaMetr_dr + 2/r )*Q,
      \
      p_dr = -( rho + p ) * ( PhiMetr_dr + alpha*Q ),
      \
      m_dr = \
        pow(r,2) * (
          4 * GV_PI * step4_A* rho
          + 1.0/2 * 1/exp_2LambdaMetr * pow(Q,2)
          + 1.0/4 * Vhat
        ),
      \
      Z_dr = \
        16*GV_PI*step4_A*(rho+p)*exp_2LambdaMetr*Omega \
        - Z * ( 4/r - (PhiMetr_dr + LambdaMetr_dr) ),
      \
      J_dr = 8.0/3*GV_PI*step4_A*(rho+p)*exp(LambdaMetr-PhiMetr)*pow(r,4)*Omega;

    dydx[1] = phiScal_dr;
    dydx[2] = Q_dr;
    dydx[3] = p_dr;
    dydx[4] = LambdaMetr_dr;
    dydx[5] = m_dr;
    dydx[6] = PhiMetr_dr;
    dydx[7] = Z_dr;
    dydx[8] = Omega_dr;
    dydx[9] = J_dr;

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

    // PhiMetr
    y0[6] = 0;

    // dOmegadr
    y0[7] = 0;

    // Omega
    y0[8] = 0;

    // J
    y0[9] = 0;

    return y0;
}

// read ODE_struct_init_modif
static double *initial_values_init_modif(void){

    double *y0 = dvector(1,COUNT_EQS);

    // phiScal
    y0[1] = phiScal_gv;

    // Q = dphiScaldr
    y0[2] = 0;

    // pressure at the center
    y0[3] = p_current;

    // Lambda - metric function
    y0[4] = 0;

    // the mass of the star
    y0[5] = 0;

    // PhiMetr
    y0[6] = 0;

    // dOmegadr
    y0[7] = 0;

    // Omega
    y0[8] = 0;

    // J
    y0[9] = 0;

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

// initializes the variables needed for odeint
// this is taking into account we are shooting only for
// PhiMetr and Omega, and thus taking into account the already shooted phiScal
static ODE_struct *ODE_struct_init_modif(void){

    ODE_struct *_ode = calloc(1, sizeof(ODE_struct));

    // set the initial values of the system
    _ode->y = initial_values_init_modif();

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

// integration function used by the shooting functions below
static void integrate_phiScal_J( ODE_struct *_ode ){

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
      &phiScal_J_foo
    );

    // let us know if the system went boom
    if(did_it_go_boom){
        printf(
          "\n integrate_phiscal says system went boom at %e \n",
          where_it_went_boom
        );

        LivePlot_phiScal_J_append(
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
        LivePlot_phiScal_J_append(
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

// integrate the ode in two steps - from start to phiScal_inf
// and from phiScal_inf to true inf with phiScal and Q set to 0
static void integrate_phiScal_J_modif( ODE_struct *_ode ){

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

    _ode->x_end = r_inf_phiscal;

    odeint(
      _ode->y, COUNT_EQS,
      _ode->x_start, _ode->x_end,
      &nok, &nbad,
      &phiScal_J_foo
    );

    // let us know if the system went boom
    if(did_it_go_boom){
        printf(
          "\n integrate_phiscal says system went boom at %e \n",
          where_it_went_boom
        );

        LivePlot_phiScal_J_append(
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
        LivePlot_phiScal_J_append(
          ODE_NAME,
          eos->model_name,
          GV_PARAMETERS_VALUES,
          R,
          COUNT_EQS
        );
    }

    _ode->x_start = r_inf_phiscal;
    _ode->x_end = r_inf;
    _ode->y[1] = 0;
    _ode->y[2] = 0;

    odeint(
      _ode->y, COUNT_EQS,
      _ode->x_start, _ode->x_end,
      &nok, &nbad,
      &phiScal_J_foo
    );

    // let us know if the system went boom
    if(did_it_go_boom){
        printf(
          "\n integrate_phiscal says system went boom at %e \n",
          where_it_went_boom
        );

        LivePlot_phiScal_J_append(
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
        LivePlot_phiScal_J_append(
          ODE_NAME,
          eos->model_name,
          GV_PARAMETERS_VALUES,
          R,
          COUNT_EQS
        );
    }

    GV_M = _ode->y[5];
    GV_J = _ode->y[9];

    odeint_point_arrs_free();

    return;
}

// the regular function which goes to regular integration process
// used for the fits integration to phiScal_inf
static void shoot_regular_execute(int n, double *v, double *f){

    // set the initial data for the ode system
    ODE_struct *_ode = ODE_struct_init();

    _ode->y[1] = v[1];
    _ode->y[6] = v[2];
    _ode->y[8] = v[3];

    integrate_phiScal_J(_ode);

    // on the left side we want scalar field to be 0
    f[1] = 0 - _ode->y[1];
    f[2] = 0 - _ode->y[6];
    f[3] = 1 - _ode->y[8];

    ODE_struct_free(&_ode);

    return;
}

// function which will use modified integrate process
// from start to phiScal and with taking into account the scalar field
// and from phiScal_inf beyond with scalar field and its derivative set to 0
// to adjust the rest shooting parameters - Omega and PhiMetr
static void shoot_regular_execute_modif(int n, double *v, double *f){

    // set the initial data for the ode system
    ODE_struct *_ode = ODE_struct_init_modif();

    _ode->y[6] = v[1];
    _ode->y[8] = v[2];

    integrate_phiScal_J_modif(_ode);

    // on the left side we want scalar field to be 0
    f[1] = 0 - _ode->y[6];
    f[2] = 1 - _ode->y[8];

    ODE_struct_free(&_ode);

    return;
}

// single shoot to provide the appropriate phiscal central value  and phiScal inf
// make another shoot to adjust the other PhiMetr inf and its central value
void single_shoot_regular_phiScal_J(void){

    // init the GV_PARAMETERS_VALUES for the scalar field
    // beta, m, lambda
    phiScal_parameters_init();

    // init the equation of state
    eos = calloc(1,sizeof(EOSmodelInfoStruct));
    eos_init(&eos);

    double v_phiScal_J = -5e-2;

    printf(
      "\n p_c = %.3e \n\n\t v_phiScal init %e \t with inf %e \n",
      p_current, v_phiScal_J, r_inf_phiscal
    );

    // lets find the interval with maximum difference
    // and the corresponding central value for the scalar field
    get_phiScal_cVal_infVal(
          GV_PARAMETERS_VALUES, p_current, eos, &v_phiScal_J,
          &r_inf_phiscal, minimal_p_power
        );

    printf(
      "\n\t got \t v_phiScal init %e \t with inf %e \n",
      v_phiScal_J, r_inf_phiscal
    );

    int newt_check = 0, newt_n = 3;

    // lets shoot for the additional indexes 6 and 8
    // in the same interval
    double *v1 = dvector(1,3);

    v1[1] = v_phiScal_J;
    v1[2] = -4e-1;
    v1[3] = 6e-1;

    r_inf = r_inf_phiscal;

    LivePlot_phiScal_J_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    newt(
        v1,
        newt_n,
        &newt_check,
        &shoot_regular_execute
    );

    printf(
      "\n\n letst try to make them better with one above %e: \t %e %e %e \n",
      r_inf, v1[1], v1[2], v1[3]
    );

    newt_n = 2;
    double *v2 = dvector(1,newt_n);

    phiScal_gv = v1[1];
    v2[1] = v1[2];
    v2[2] = v1[3];
    r_inf = 1e4;

    newt(
        v2,
        newt_n,
        &newt_check,
        &shoot_regular_execute_modif
    );

    printf(
      "\n\n New %e: \t %e %e \n",
      r_inf, v2[1], v2[2]
    );

    free(v1);
    free(v2);
    eos_free(&eos);
    free(GV_PARAMETERS_VALUES);

    return;
}

// single shoot to provide the appropriate phiscal central value  and phiScal inf
// make another shoot to adjust the other PhiMetr inf and its central value
// it will increase the infinity until no significant change occurs in
// PhiMetr_c and Omega_c
void single_shoot_regular_phiScal_J_iterate_inf(void){

    // init the GV_PARAMETERS_VALUES for the scalar field
    // beta, m, lambda
    phiScal_parameters_init();

    // init the equation of state
    eos = calloc(1,sizeof(EOSmodelInfoStruct));
    eos_init(&eos);

    double v_phiScal_J = -5e-2;

    printf(
      "\n p_c = %.3e \n\n\t v_phiScal init %e \t with inf %e \n",
      p_current, v_phiScal_J, r_inf_phiscal
    );

    // lets find the interval with maximum difference
    // and the corresponding central value for the scalar field
    get_phiScal_cVal_infVal(
          GV_PARAMETERS_VALUES, p_current, eos, &v_phiScal_J,
          &r_inf_phiscal, minimal_p_power
        );

    printf(
      "\n\t got \t v_phiScal init %e \t with inf %e \n",
      v_phiScal_J, r_inf_phiscal
    );

    int newt_check = 0, newt_n = 3;

    // lets shoot for the additional indexes 6 and 8
    // in the same interval
    double *v1 = dvector(1,3);

    v1[1] = v_phiScal_J;
    v1[2] = -4e-1;
    v1[3] = 6e-1;

    r_inf = r_inf_phiscal;

    newt(
        v1,
        newt_n,
        &newt_check,
        &shoot_regular_execute
    );

    printf(
      "\n\n letst try to make them better with one above %e: \t %e %e %e \n",
      r_inf, v1[1], v1[2], v1[3]
    );

    newt_n = 2;
    double \
      *v2 = dvector(1,newt_n),
      PhiMetr_c = v1[2],
      Omega_c = v1[3];

    phiScal_gv = v1[1];
    v2[1] = PhiMetr_c;
    v2[2] = Omega_c;
    r_inf = 1e4;

    int count = 0;
    while(++count){

        newt(
            v2,
            newt_n,
            &newt_check,
            &shoot_regular_execute_modif
        );

        printf("\n\t iteration %d: %e %e for inf %e \n", count, v2[1], v2[2], r_inf);

        if(get_power(PhiMetr_c - v2[1]) >= -6 || get_power(Omega_c - v2[2]) >= -6 ){
            PhiMetr_c = v2[1];
            Omega_c = v2[2];
            r_inf += 0.2*r_inf;
        }else{
            printf("\n no sig diff \n");
            break;
        }
    }

    printf(
      "\n\n New %e: \t %e %e \n",
      r_inf, v2[1], v2[2]
    );

    LivePlot_phiScal_J_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    ODE_struct *_ode = ODE_struct_init_modif();
    _ode->y[6] = v2[1];
    _ode->y[8] = v2[2];

    integrate_phiScal_J_modif(_ode);

    // stop for further investigation if the boundary condition for the
    // scalar field is not satisfied
    // should be connected to the solver_newt, but I prefer to set it here
    // as the power of -10 by hand
    if(get_power(0 - _ode->y[1]) > -9 || get_power(0 - _ode->y[6]) > -9 || get_power(1 - _ode->y[8]) > -9 ){
        printf("\n scalar field boundary condition not met, terminating...");
        exit(123);
    }

    ODE_struct_free(&_ode);
    free(v1);
    free(v2);
    eos_free(&eos);
    free(GV_PARAMETERS_VALUES);

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

    }else if( get_power(*v_phiScal + tiny) > -2 ){
        printf(
          "\n\t\t changing the initial guess from %e to %e \n",
          *tmp_phiScal,
          *v_phiScal
        );

        *tmp_phiScal = *v_phiScal;
    }

    return;
}

// iterate over the pressures
// single shoot to provide the appropriate phiscal central value  and phiScal inf
// make another shoot to adjust the other PhiMetr inf and its central value
// it will increase the infinity until no significant change occurs in
// PhiMetr_c and Omega_c
void single_shoot_regular_phiScal_J_iterate_inf_iterpres(void){

    // init the GV_PARAMETERS_VALUES for the scalar field
    // beta, m, lambda
    phiScal_parameters_init();

    // init the equation of state
    eos = calloc(1,sizeof(EOSmodelInfoStruct));
    eos_init(&eos);

    double \
      // initial guessed values
      v_phiScal_J = -5e-2,
      v_phiScal_J_tmp = v_phiScal_J,
      PhiMetr_c = -4e-1,
      Omega_c = 6e-1,
      // central pressure values
      pressure_start = P_START, pressure_end = P_END;

    p_current = pressure_start;

    // resetting the result file
    ResultFile_phiScal_J_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

    // start going over the pressures starting with
    // P_START until reaching the P_END
    // and changing it in a specific way in change_current_pressure()
    int iterate = 1;
    while(iterate){

        // reset the file for live plotting
        LivePlot_phiScal_J_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

        // set initial guess for the infinity for the scalar field
        r_inf_phiscal = R_INF_PHISCAL;

        printf("\n p_c = %.3e \n",p_current);
        v_phiScal_J = v_phiScal_J_tmp;
        printf(
          "\n\t scalar field shoot with initial"
          "v_phiScal = %e and phiScal_inf = %e \n",
          v_phiScal_J, r_inf_phiscal
        );

        // lets find the interval with maximum difference
        // and the corresponding central value for the scalar field
        get_phiScal_cVal_infVal(
          GV_PARAMETERS_VALUES, p_current, eos, &v_phiScal_J,
          &r_inf_phiscal, minimal_p_power
        );

        printf(
          "\n\t setting v_phiScal = %e and phiScal_inf = %e \n",
          v_phiScal_J, r_inf_phiscal
        );

        int newt_check = 0;
        // lets shoot for the additional indexes 6 and 8
        // in the same interval in which no significant difference
        // occurs in the initial scalar field
        double *v1 = dvector(1,3);
        v1[1] = v_phiScal_J;
        v1[2] = PhiMetr_c;
        v1[3] = Omega_c;

        r_inf = r_inf_phiscal;

        printf(
            "\n\t fine shoot for inf %e \n"
            "\n\t\t ( %e %e %e ) ---> ",
            r_inf_phiscal, v1[1], v1[2], v1[3]
        );

        newt(
            v1,
            3,
            &newt_check,
            &shoot_regular_execute
        );

        printf(
            " ( %e %e %e ); R = %e \n",
            v1[1], v1[2], v1[3], R
        );

        v_phiScal_J = v1[1];
        PhiMetr_c = v1[2];
        Omega_c = v1[3];
        double *v2 = dvector(1,2), inftmp = R_INF;

        // set the guesses in the second integrator who
        // changes the scalar field and its derivative to 0
        // when it reaches the r_inf_phiscal
        // and it shoots for the PhiMetr and Omega_c with increasing r_inf
        // until no significant difference occurs
        phiScal_gv = v1[1];
        v2[1] = PhiMetr_c;
        v2[2] = Omega_c;
        r_inf = R_INF;

        int count = 0;
        while(++count){

            newt(
                v2,
                2,
                &newt_check,
                &shoot_regular_execute_modif
            );

            printf(
                "\n\t v_%d ---> for inf %e ( R = %e )"
                "\n\t\t %e ---> %e  \n"
                "\n\t\t %e ---> %e  \n",
                count, r_inf, R,
                PhiMetr_c, v2[1],
                Omega_c, v2[2]
            );

            // check whether there is any significant difference
            // in the now and previous shooted central values in
            // the PhiMetr_c and Omega_c
            if((
                get_power(PhiMetr_c - v2[1]) > -10 || \
                get_power(Omega_c - v2[2]) > -10
            )){
                PhiMetr_c = v2[1];
                Omega_c = v2[2];
                inftmp = r_inf;
                r_inf += 0.2*r_inf;
            }else{
                printf(
                    "\n\t No sig diff in with \n\t\t %e %e and inf %e \n",
                    v2[1], v2[2], r_inf
                );
                break;
            }
        }

        // after all this effort append the results
        // the mass and inertial moment are saved at
        // integrate_phiScal_J_modif()
        ResultFile_phiScal_J_append(
          ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES,
          p_current, v_phiScal_J, GV_M,
          AR,
          rho_c, GV_J
        );

        free(v1);
        free(v2);

        // check if there is some scalarsization or not
        // to carry the guessed value or not
        check_achieved_phiScal(&v_phiScal_J, &v_phiScal_J_tmp);

        // increment the central pressure power wise
        change_current_pressure();
        // if the current central pressure is greater the greatest desired one
        // just end the loop
        if(p_current > pressure_end){
            iterate = 0;
        }else{
            iterate = 1;
        }

        // some integrations at the end to take data
        // and make profiles
        // but due to the stiffnes it did not worked well
        // thus commeting out and maybe getting rid of it
        // after some commit
        //ODE_struct *_ode2 = ODE_struct_init();
        //_ode2->y[1] = v_phiScal_J;
        //_ode2->y[6] = v2[1];
        //_ode2->y[8] = v2[2];
        //_ode2->x_end = R - 1e-6*R;

        //double tmpAR = AR, small_step = 1e-7;
        //integrate_phiScal_J(_ode2);

        // lets try to find the radii by hand
        //while(get_power(_ode2->y[3]) > minimal_p_power){
            //_ode2->x_start = _ode2->x_end;
            //_ode2->x_end += small_step;
            //integrate_phiScal_J(_ode2);
            //printf("\n %.7e --- %.7e \n",_ode2->x_end, _ode2->y[3]);
        //}

        //LivePlot_phiScal_J_open(ODE_NAME, eos->model_name, GV_PARAMETERS_VALUES);

        //ODE_struct *_ode = ODE_struct_init();
        //_ode->y[1] = v_phiScal_J;
        //_ode->y[6] = v2[1];
        //_ode->y[8] = v2[2];
        //_ode->x_end = r_inf_phiscal;

        //printf(
            //"\n\t\t final integr with \n\t inf %e \t %e %e %e \n",
            //_ode->x_end, _ode->y[1], _ode->y[6], _ode->y[8]
        //);

        //integrate_phiScal_J(_ode);

        //if(get_power(0 - _ode->y[1]) > -10){
            //printf("\n scalar field boundary condition not met, terminating...");
            //exit(123);
        //}

        //_ode->y[1] = 0;
        //_ode->y[2] = 0;
        //_ode->x_start = r_inf_phiscal;
        //_ode->x_end = r_inf;

        //printf(
            //"\n\t\t final integr with \n\t inf %e \t %e %e %e \n",
            //_ode->x_end, _ode->y[1], _ode->y[6], _ode->y[8]
        //);

        //integrate_phiScal_J(_ode);

        // stop for further investigation if the boundary condition for the
        // scalar field is not satisfied
        // should be connected to the solver_newt, but I prefer to set it here
        // as the power of -10 by hand
        //if( get_power(0 - _ode->y[6]) > -10 || get_power(1 - _ode->y[8]) > -10 ){
            //printf("\n other two boundary condition not met, terminating...");
            //exit(123);
        //}

        //ODE_struct_free(&_ode);
        //ODE_struct_free(&_ode2);
        //free(v1);
        //free(v2);

        //check_achieved_phiScal(&v_phiScal_J, &v_phiScal_J_tmp);

        // increment the central pressure power wise
        //change_current_pressure();
        // if the current central pressure is greater the greatest desired one
        // just end the loop
        //if(p_current > pressure_end){
            //iterate = 0;
        //}else{
            //iterate = 1;
        //}
    }

    eos_free(&eos);
    free(GV_PARAMETERS_VALUES);

    return;
}
