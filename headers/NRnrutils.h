#ifndef _NRnrutils_HEADERS_H
    #define _NRnrutils_HEADERS_H

    #include <stddef.h>
    #define NR_END 1
    #define FREE_ARG char*

    extern double sqrarg;
    #define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
    extern double dsqrarg;
    #define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
    extern double dmaxarg1,dmaxarg2;
    #define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
    extern double dminarg1,dminarg2;
    #define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
    extern double maxarg1,maxarg2;
    #define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
    extern double minarg1,minarg2;
    #define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
    extern long lmaxarg1,lmaxarg2;
    #define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))
    extern long lminarg1,lminarg2;
    #define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2))
    extern int imaxarg1,imaxarg2;
    #define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))
    extern int iminarg1,iminarg2;
    #define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
    #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

    extern void nrerror(const char error_text[]);
    extern float *vector(long nl, long nh);
    extern int *ivector(long nl, long nh);
    extern unsigned char *cvector(long nl, long nh);
    extern unsigned long *lvector(long nl, long nh);
    extern double *dvector(long nl, long nh);
    extern float **matrix( \
        unsigned long nrl, \
        unsigned long nrh, \
        unsigned long ncl, \
        unsigned long nch \
    );
    extern double **dmatrix(long nrl, long nrh, long ncl, long nch);
    extern int **imatrix(long nrl, long nrh, long ncl, long nch);
    extern float **submatrix( \
        float **a, \
        unsigned long oldrl, \
        unsigned long oldrh, \
        unsigned long oldcl, \
        unsigned long oldch, \
        unsigned long newrl, \
        unsigned long newcl \
    );
    extern float **convert_matrix(\
        float *a, \
        unsigned long nrl, \
        unsigned long nrh, \
        unsigned long ncl, \
        unsigned long nch \
    );
    extern float ***f3tensor(\
        unsigned long nrl, \
        unsigned long nrh, \
        unsigned long ncl, \
        unsigned long nch, \
        unsigned long ndl, \
        unsigned long ndh \
    );
    extern void free_vector(float *v, long nl, long nh);
    extern void free_ivector(int *v, long nl, long nh);
    extern void free_cvector(unsigned char *v, long nl, long nh);
    extern void free_lvector(unsigned long *v, long nl, long nh);
    extern void free_dvector(double *v, long nl, long nh);
    extern void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
    extern void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
    extern void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
    extern void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
    extern void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
    extern void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);

#endif

