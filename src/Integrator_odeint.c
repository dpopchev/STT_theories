#include "ExternalHeaders.h"

// odeint defines
#define EPS 1e-12
#define H1 1e-30
#define HMIN 0.0
#define MAXSTP 10000
#define TINY 1.0e-30

// rkqs defines
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

// odeint determines the step size by monitoring
// local truncation error by scaling several scaling techniques are available
// thus choose one of the below
#define ODEINT_SCALING_METHOD 0

// rkqs is advancing with small steps
// choose the method to evaluate the new value for x
#define RKQS_STEP_METHOD 0

static const char *ODEINT_SCALING_METHOD_DESCRIPTION[] = \
    { \
        "[0], \
         the default, \
         yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;", \
         \
        "[1], \
         alternative scaling from the NR, \
         yscal[i]=fabs(eps*dydx[i]*h)+TINY;", \
        \
        "[2], \
         stiff integrator recommends scaling to constant value, try 1, \
         yscal[i]=FMAX(fabs(eps*dydx[i]*h), 1e-11);", \
        \
        "[3], \
         stiff integrator recommends scaling to constant value, try 2, \
         yscal[i]=FMAX(fabs(y[i]), 1e-15);", \
        \
        "[4], \
         stiff integrator recommends scaling to constant value, try 3, \
         yscal[i]=FMAX(fabs(y[i])+fabs(dydx[i]*h)+TINY, 1e-11);" \
    };

static const char *RKQS_STEP_METHOD_DESCRIPTION[] = \
    { \
        "[0], \
         the default method, \
         xnew = (*x) + h;", \
        \
        "[1], \
         borrowed from the numerical diff methods, \
         volatile double hh = *x + h; \
         hh -= *x; \
         xnew = (*x) + hh;"\
    };

void odeint_info_print_stdout(void){

    printf(
        "\n Integrator info: \n"
        "\n\t General info: \n"
        "\t\t name: %s \n"
        "\n\t Scaling method: \n"
        "\t\t %s \n"
        "\n\t RKQS step method: \n"
        "\t\t %s \n"
        "\n\t Parameters: \n"
        "\t\t initial step %.3e \n"
        "\t\t minimal step %.3e \n"
        "\t\t desired accuracy %.3e \n"
        "\t\t safety parameter is %.3e \n"
        "\t\t pgrow parameter is %.3e \n"
        "\t\t pshrink parameter is %.3e \n"
        "\t\t errcon parameter is %.3e \n"
        "\t\t maxstp parameter is %d \n"
        "\t\t tiny parameter is %.3e \n"
        ,"odeint",
        ODEINT_SCALING_METHOD_DESCRIPTION[ODEINT_SCALING_METHOD],
        RKQS_STEP_METHOD_DESCRIPTION[RKQS_STEP_METHOD],
        H1, HMIN, EPS, SAFETY, PGROW, PSHRNK, ERRCON, MAXSTP, TINY
    );

    return;
}

void odeint_info_print_ResultFile(FILE *fp){

    fprintf(
        fp,
        "\n Integrator info: \n"
        "\n\t General info: \n"
        "\t\t name: %s \n"
        "\n\t Scaling method: \n"
        "\t\t %s \n"
        "\n\t RKQS step method: \n"
        "\t\t %s \n"
        "\n\t Parameters: \n"
        "\t\t initial step %.3e \n"
        "\t\t minimal step %.3e \n"
        "\t\t desired accuracy %.3e \n"
        "\t\t safety parameter is %.3e \n"
        "\t\t pgrow parameter is %.3e \n"
        "\t\t pshrink parameter is %.3e \n"
        "\t\t errcon parameter is %.3e \n"
        "\t\t maxstp parameter is %d \n"
        "\t\t tiny parameter is %.3e \n"
        ,"odeint",
        ODEINT_SCALING_METHOD_DESCRIPTION[ODEINT_SCALING_METHOD],
        RKQS_STEP_METHOD_DESCRIPTION[RKQS_STEP_METHOD],
        H1, HMIN, EPS, SAFETY, PGROW, PSHRNK, ERRCON, MAXSTP, TINY
    );

    fclose(fp);

    return;
}

// odeint global variables
int kmax,kount;
double *xp,**yp,dxsav;

// rkck code
static void rkck(\
    double y[], double dydx[], int n, double x, double h, double yout[],
    double yerr[], void (*derivs)(double, double [], double [])\
){
    int i;
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
        b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
        b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
        b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
        b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
        c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
        dc5 = -277.00/14336.0;
    double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
        dc4=c4-13525.0/55296.0,dc6=c6-0.25;
    double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

    ak2=dvector(1,n);
    ak3=dvector(1,n);
    ak4=dvector(1,n);
    ak5=dvector(1,n);
    ak6=dvector(1,n);
    ytemp=dvector(1,n);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
    (*derivs)(x+a3*h,ytemp,ak3);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    (*derivs)(x+a4*h,ytemp,ak4);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    (*derivs)(x+a5*h,ytemp,ak5);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    (*derivs)(x+a6*h,ytemp,ak6);
    for (i=1;i<=n;i++)
        yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=1;i<=n;i++)
        yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    free_dvector(ytemp,1,n);
    free_dvector(ak6,1,n);
    free_dvector(ak5,1,n);
    free_dvector(ak4,1,n);
    free_dvector(ak3,1,n);
    free_dvector(ak2,1,n);
}

// rkqs code
static void rkqs(\
    double y[], double dydx[], int n, double *x, double htry, double eps,
    double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double []) \
){
    //void rkck(
        //double y[], double dydx[], int n, double x, double h,
        //double yout[], double yerr[],
        //void (*derivs)(double, double [], double [])
    //);
    int i;
    double errmax,h,htemp,xnew,*yerr,*ytemp;

    yerr=dvector(1,n);
    ytemp=dvector(1,n);
    h=htry;
    for (;;) {
        rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
        errmax=0.0;
        for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
        errmax /= eps;
        if (errmax <= 1.0) break;
        htemp=SAFETY*h*pow(errmax,PSHRNK);
        h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));

        switch(RKQS_STEP_METHOD){
            case 0:
                xnew = (*x) + h;
                break;
            case 1: {
                volatile double hh = *x + h;
                hh -= *x;
                xnew = (*x) + hh;
                break; }
            default:
                printf(\
                    "\n rkqs step method %d UNKNOWN, line 213", \
                    RKQS_STEP_METHOD\
                );
                    exit(213);
                    break;
        }
        xnew=(*x)+h;
        if (xnew == *x) nrerror("stepsize underflow in rkqs");
    }
    if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
    else *hnext=5.0*h;
    *x += (*hdid=h);
    for (i=1;i<=n;i++) y[i]=ytemp[i];
    free_dvector(ytemp,1,n);
    free_dvector(yerr,1,n);
}

// odeint code
void odeint(
    double ystart[], int nvar, double x1, double x2, int *nok, int *nbad,
    void (*derivs)(double, double [], double [])
){
    int nstp,i;
    double xsav = 0, x,hnext,hdid,h,eps=EPS, h1=H1, hmin=HMIN;
    double *yscal,*y,*dydx;

    yscal=dvector(1,nvar);
    y=dvector(1,nvar);
    dydx=dvector(1,nvar);
    x=x1;
    h=SIGN(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1;i<=nvar;i++) y[i]=ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;
    for (nstp=1;nstp<=MAXSTP;nstp++) {
        (*derivs)(x,y,dydx);

        for (i=1;i<=nvar;i++){
            switch(ODEINT_SCALING_METHOD){
                case 0:
                    yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
                    break;
                case 1:
                    yscal[i]=fabs(eps*dydx[i]*h)+TINY;
                    break;
                case 2:
                    yscal[i]=FMAX(fabs(eps*dydx[i]*h), 1e-11);
                    break;
                case 3:
                    yscal[i]=FMAX(fabs(y[i]), 1e-15);
                    break;
                case 4:
                    yscal[i]=FMAX(fabs(y[i])+fabs(dydx[i]*h)+TINY, 1e-11);
                    break;
                default:
                    printf(\
                        "\n odeint.c odeint scaling method %d UNKNOWN, line 269", \
                        ODEINT_SCALING_METHOD\
                    );
                    exit(269);
                    break;
            }
        }

        if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
            xp[++kount]=x;
            for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
            xsav=x;
        }
        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
        (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
        if (hdid == h) ++(*nok); else ++(*nbad);
        if ((x-x2)*(x2-x1) >= 0.0) {
            for (i=1;i<=nvar;i++) ystart[i]=y[i];
            if (kmax) {
                xp[++kount]=x;
                for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
            }
            free_dvector(dydx,1,nvar);
            free_dvector(y,1,nvar);
            free_dvector(yscal,1,nvar);
            return;
        }
        if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
        h=hnext;
    }
    nrerror("Too many steps in routine odeint");
}

#undef EPS
#undef H1
#undef HMIN
#undef MAXSTP
#undef TINY
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef ODEINT_SCALING_METHOD
#undef RKQS_STEP_METHOD
