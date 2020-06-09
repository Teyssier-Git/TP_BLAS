#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

#define abs(a) ((a) > 0 ? (a) : (-a))

CBLAS_INDEX mnblas_isamax(const int n, const float *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    CBLAS_INDEX currentMax = 0;
    float maxValue = abs(x[0]);

    for(int i=1; i<n; i += incx){
        float tried = abs(x[i]);
        if(tried > maxValue){
            currentMax = i;
            maxValue = tried;
        }
    }

    return currentMax+1;

}

CBLAS_INDEX mnblas_idamax(const int n, const double *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    CBLAS_INDEX currentMax = 0;
    double maxValue = abs(x[0]);

    for(int i=1; i<n; i += incx){
        double tried = abs(x[i]);
        if(tried > maxValue){
            currentMax = i;
            maxValue = tried;
        }
    }

    return currentMax+1;

}

CBLAS_INDEX mnblas_icamax(const int n, const void *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    complexe_float_t *X = (complexe_float_t *)x;

    CBLAS_INDEX currentMax = 0;
    float maxValue = abs(X[0].real) + abs(X[0].imaginary);

    for(int i=1; i<n; i += incx){
        float tried = abs(X[i].real) + abs(X[i].imaginary);
        if(tried > maxValue){
            currentMax = i;
            maxValue = tried;
        }
    }

    return currentMax+1;

}

CBLAS_INDEX mnblas_izamax(const int n, const void *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    complexe_double_t *X = (complexe_double_t *)x;

    CBLAS_INDEX currentMax = 0;
    double maxValue = abs(X[0].real) + abs(X[0].imaginary);

    for(int i=1; i<n; i += incx){
        double tried = abs(X[i].real) + abs(X[i].imaginary);
        if(tried > maxValue){
            currentMax = i;
            maxValue = tried;
        }
    }

    return currentMax+1;

}
