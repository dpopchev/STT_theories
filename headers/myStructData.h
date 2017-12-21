#ifndef _myStruct_DATA_H_
    #define _myStruct_DATA_H_

    // structure to describe the ODE system
    // and the integrator to be used
    typedef struct ODEsystemStruct{
        int \
            // ODE system number - description and functions
            index, \
            // total number of ODE equations
            eqs_count, \
            // amount of points we want to save, if any, see below
            points_count, \
            // count OK steps done by the integrator
            nok, \
            // count BAD steps done by the integrator
            nbad, \
            // scaling method (odeint), see Integrator_odeint.c for more details
            // default is 0
            odeint_scaling_method, \
            // method to evaluate the new value for x in rkqs
            // see odeint.c for more details
            // default is 0
            rkqs_step_method;

        double \
            // the system itself
            *y, \
            // the independent variable start
            x_initial, \
            // the independent variable
            x_final, \
            // place to save x points, size is defined by points_count
            *points_x, \
            // place to save corresponding y points, size is defined by points_count
            **points_y, \
            // desired accuracy for the integrator (odeint)
            eps, \
            // initial step size for the integrator (odeint)
            h1, \
            // minimum allowed step size, can be zero (odeint)
            hmin;

        char \
            // name with maybe simple description, see GlobalVars.h
            name_system[512], \
            // list variables index, see GlobalVars.h
            name_vars[512], \
            // odeint scaling method used explicit written
            odeint_scaling_method_description[512], \
            // rkqs advancing method in odeint used explicit written
            rkqs_step_method_description[512];

        void (*foo)(double, double *, double *);

    } ODEsystemStruct;

#endif
