#include "ExternalHeaders.h"

#define DEBUGGING 1

void dvector_copy( double *copy_from, double *copy_to, int size ){

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
                "%s %s copy_from[%d] = %.3e copy_to[%d] = %.3e \n", \
                identation, function_path, \
                i, copy_from[i], \
                i, copy_to[i] \
            );
        }

        copy_to[i] = copy_from[i];
    }

    return;
}

#undef DEBUGGING
