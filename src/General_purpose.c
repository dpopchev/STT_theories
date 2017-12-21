#include "ExternalHeaders.h"

#define DEBUGGING 0

void dvector_copy( double *from, double *to, int size ){

    const char \
        function_path[] = "General_purpose.c dvector_copy", \
        identation[] = "\n";

    if(DEBUGGING){
        printf(\
            "%s %s starting to copy dvector with size %d \n", \
            identation, function_path, size \
        );
    }

    for(int i=1; i <= size; i++){

        if(DEBUGGING){
            printf(\
                "%s %s from[%d] = %.3e to[%d] = %.3e \n", \
                identation, function_path, \
                i, from[i], \
                i, to[i] \
            );
        }

        to[i] = from[i];
    }

    return;
}

#undef DEBUGGING
