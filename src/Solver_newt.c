#include "ExternalHeaders.h"

// newt variable for maximum number of iterations
#define MAXITS 200

// newt sets convergence criterion on function values
#define TOLF 1.0e-4

// newt sets the criterion for deciding whether spurious convergence to a minimum
// of NR_fmin has occured
#define TOLMIN 1.0e-6

// newt is the convergence criterion on \delta x, such that it is not meaningful to
// make more corrections to the root
#define TOLX_newt 1.0e-7

// newt maximum scaled step length allowed in the line searches
#define STPMX 100.0

// lnsrch ensures sufficient decrease in function value
#define ALF 1.0e-14

// lnsrch convergence criterion on \Delta x
#define TOLX_lnsrch 1.0e-12

// machine epsilon for the numerical Jacobian fdjac
#define EPS_fdjac 1.0e-4

// ludcmp smallest number
#define TINY 1e-30

static int nn;
static double *fvec;
static void (*nrfuncv)(int n, double v[], double f[]);
#define FREERETURN {free_dvector(fvec,1,n);free_dvector(xold,1,n);\
    free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);\
    free_ivector(indx,1,n);return;}

static void ludcmp(double **a, int n, int *indx, double *d){

    int i,imax=0,j,k;
    double big,dum,sum,temp;
    double *vv;

    vv=dvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,n);
}

static void lubksb(double **a, int n, int *indx, double b[]){
    int i,ii=0,ip,j;
    double sum;

    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}

static void lnsrch(
    int n, double xold[], double fold, double g[], double p[],
    double x[], double *f, double stpmax, int *check,
    double (*func)(double [])
){

    int i;
    double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
        test,tmplam;

    *check=0; f2 = 0; alam2 = 0;

    for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
    sum=sqrt(sum);
    if (sum > stpmax)
        for (i=1;i<=n;i++) p[i] *= stpmax/sum;
    for (slope=0.0,i=1;i<=n;i++)
        slope += g[i]*p[i];
    if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
    test=0.0;
    for (i=1;i<=n;i++) {
        temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
        if (temp > test) test=temp;
    }
    alamin=TOLX_lnsrch/test;
    alam=1.0;
    for (;;) {
        for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
        *f=(*func)(x);
        if (alam < alamin) {
            for (i=1;i<=n;i++) x[i]=xold[i];
            *check=1;
            return;
        } else if (*f <= fold+ALF*alam*slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope/(2.0*(*f-fold-slope));
            else {
                rhs1 = *f-fold-alam*slope;
                rhs2=f2-fold-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else {
                    disc=b*b-3.0*a*slope;
                    if (disc < 0.0) tmplam=0.5*alam;
                    else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
                    else tmplam=-slope/(b+sqrt(disc));
                }
                if (tmplam > 0.5*alam)
                    tmplam=0.5*alam;
            }
        }
        alam2=alam;
        f2 = *f;
        alam=FMAX(tmplam,0.1*alam);
    }
}

// numerical Jacobian using forward difference approximation method
static void fdjac(int n, double x[], double _fvec[], double **df,
    void (*vecfunc)(int, double [], double [])
){
    int i,j;
    double h,temp,*f;

    f=dvector(1,n);
    for (j=1;j<=n;j++) {
        temp=x[j];
        h=EPS_fdjac*fabs(temp);
        if (h == 0.0) h=EPS_fdjac;
        x[j]=temp+h;
        h=x[j]-temp;
        (*vecfunc)(n,x,f);
        x[j]=temp;
        for (i=1;i<=n;i++) df[i][j]=(f[i]-_fvec[i])/h;
    }
    free_dvector(f,1,n);
}

// returns f = 0.5 \vec F * \vec F at x
static double NR_fmin(double x[]){
    int i;
    double sum;

    (*nrfuncv)(nn,x,fvec);
    for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
    return 0.5*sum;
}

void newt_info_print_stdout(void){

    printf(
        "\n Solver info: \n"
        "\n\t General info \n"
        "\t\t name %s \n"
        "\n\t %s parameters info \n"
        "\t\t max iterations %d \n"
        "\t\t function value convergence criteria %.3e \n"
        "\t\t spurious convergence to minimum criteria %.3e \n"
        "\t\t criteria for meaningful corrections to root %.3e \n"
        "\t\t maximum sclaed step length in line searches %.3e \n"
        "\n\t lnsrch parameters \n"
        "\t\t sufficient decrease in function value %.3e \n"
        "\t\t convergence criterion on Delta x %.3e \n"
        "\n\t numerical Jacobian fdjac \n"
        "\t\t differentiation method: forward difference method \n"
        "\t\t machine epsilon %.3e \n"
        "\n\t ludcmp parameters \n"
        "\t\t tiny parameter %.3e \n"
        , "newt","newt",
        MAXITS, TOLF, TOLMIN, TOLX_newt, STPMX, ALF, TOLX_lnsrch,
        EPS_fdjac, TINY
    );

    return;
}

void newt_info_print_ResultFile(FILE *fp){

    fprintf(
        fp,
        "\n Solver info: \n"
        "\n\t General info \n"
        "\t\t name %s \n"
        "\n\t %s parameters info \n"
        "\t\t max iterations %d \n"
        "\t\t function value convergence criteria %.3e \n"
        "\t\t spurious convergence to minimum criteria %.3e \n"
        "\t\t criteria for meaningful corrections to root %.3e \n"
        "\t\t maximum sclaed step length in line searches %.3e \n"
        "\n\t lnsrch parameters \n"
        "\t\t sufficient decrease in function value %.3e \n"
        "\t\t convergence criterion on Delta x %.3e \n"
        "\n\t numerical Jacobian fdjac \n"
        "\t\t differentiation method: forward difference method \n"
        "\t\t machine epsilon %.3e \n"
        "\n\t ludcmp parameters \n"
        "\t\t tiny parameter %.3e \n"
        , "newt","newt",
        MAXITS, TOLF, TOLMIN, TOLX_newt, STPMX, ALF, TOLX_lnsrch,
        EPS_fdjac, TINY
    );

    fclose(fp);

    return;
}

void newt(double x[], int n, int *check, void (*vecfunc)(int, double [], double [])){
    //void fdjac(int n, double x[], double fvec[], double **df,
        //void (*vecfunc)(int, double [], double []));
    //double fmin(double x[]);
    //void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
         //double *f, double stpmax, int *check, double (*func)(double []));
    //void lubksb(double **a, int n, int *indx, double b[]);
    //void ludcmp(double **a, int n, int *indx, double *d);

    int i,its,j,*indx;
    double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

    indx=ivector(1,n);
    fjac=dmatrix(1,n,1,n);
    g=dvector(1,n);
    p=dvector(1,n);
    xold=dvector(1,n);
    fvec=dvector(1,n);
    nn=n;
    nrfuncv=vecfunc;
    f=NR_fmin(x);
    test=0.0;
    for (i=1;i<=n;i++)
        if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < 0.01*TOLF) {
        *check=0;
        FREERETURN
    }
    for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
    stpmax=STPMX*FMAX(sqrt(sum),(double)n);
    for (its=1;its<=MAXITS;its++) {
        fdjac(n,x,fvec,fjac,vecfunc);
        for (i=1;i<=n;i++) {
            for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
            g[i]=sum;
        }
        for (i=1;i<=n;i++) xold[i]=x[i];
        fold=f;
        for (i=1;i<=n;i++) p[i] = -fvec[i];
        ludcmp(fjac,n,indx,&d);
        lubksb(fjac,n,indx,p);
        lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,NR_fmin);
        test=0.0;
        for (i=1;i<=n;i++)
            if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
        if (test < TOLF) {
            *check=0;
            FREERETURN
        }
        if (*check) {
            test=0.0;
            den=FMAX(f,0.5*n);
            for (i=1;i<=n;i++) {
                temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
                if (temp > test) test=temp;
            }
            *check=(test < TOLMIN ? 1 : 0);
            FREERETURN
        }
        test=0.0;
        for (i=1;i<=n;i++) {
            temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
            if (temp > test) test=temp;
        }
        if (test < TOLX_newt) FREERETURN
    }
    nrerror("MAXITS exceeded in newt");
}

#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX_newt
#undef STPMX
#undef FREERETURN
#undef EPS_fdjac
#undef ALF
#undef TOLX
#undef TINY
