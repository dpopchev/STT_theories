#ifndef _myStruct_DATA_H_
    #define _myStruct_DATA_H_

    // structure to describe the ODE system
    typedef struct ODEsystemStruct{
        int \
            // ODE system number - description and functions
            index,
            // total number of ODE equations
            eqs_count,
            // amount of points we want to save, if any, see below
            points_count,
            // number of OK and BAD steps done by the integrator
            nok, nbad;

        double \
            // the system itself
            *y, \
            // the independent variable
            x, \
            // if we want to save separate points we will use this
            // see Numerical Recepies
            *points_x, **points_y;

        char \
            // name with maybe simple description, see GlobalVars.h
            name_system[512], \
            // list variables index, see GlobalVars.h
            name_vars[512];

        void (*foo)(double, double *, double *);

    } ODEsystemStruct;

#endif
