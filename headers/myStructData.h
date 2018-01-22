#ifndef _myStruct_DATA_H_
    #define _myStruct_DATA_H_

    // structure to describe the ODE system
    // and the integrator to be used
    typedef struct ODEsystemStruct{
        int \
            // ODE system index
            index, \
            // total number of ODE equations
            eqs_count, \
            // amount of points we want to record, if any
            points_count, \
            // count OK steps done by the integrator
            nok, \
            // count BAD steps done by the integrator
            nbad, \
            // scaling method, see Integrator_odeint.c for more details
            // default is 0
            odeint_scaling_method, \
            // method to evaluate the new value for x in rkqs
            // see Integrator_odeint.c for more details
            // default is 0
            rkqs_step_method, \
            // amount of free parameters for the system
            free_parmeters_count_all, \
            // how many values of each parameter we are interested
            *free_parameters_count_each, \
            // index for the y whos initial value will be changed
            index_of_y_to_change;

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
            hmin, \
            // pointer to include the free parameters of the system
            *free_parmeters_values, \
            // all possible values of parameters we are interested per parameter
            **free_parmeters_values_all, \
            // interval in which we will change the initial value of y
            // with index index_of_y_to_change
            initial_y_start, initial_y_end, initial_y_current, initial_y_step;

        char \
            // name with simple description, see GlobalVars.h
            name_system[512], \
            // list variables index, see GlobalVars.h
            **name_vars, \
            // odeint scaling method used explicit written
            odeint_scaling_method_description[512], \
            // rkqs advancing method in odeint used explicit written
            rkqs_step_method_description[512], \
            // names and some description of the free parameters of the system
            **free_parmeters_names;

        void (*foo)(double, double *, double *);

    } ODEsystemStruct;

#endif
