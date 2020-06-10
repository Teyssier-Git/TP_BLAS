#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

#define abs(a) ((a) > 0 ? (a) : (-a))

CBLAS_INDEX mnblas_isamin(const int n, const float *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    CBLAS_INDEX currentMin = 0;
    float minValue = abs(x[0]);

    for(int i=1; i<n; i += incx){
        float tried = abs(x[i]);
        if(tried < minValue){
            currentMin = i;
            minValue = tried;
        }
    }

    return currentMin+1;

}

CBLAS_INDEX mnblas_idamin(const int n, const double *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    CBLAS_INDEX currentMin = 0;
    double minValue = abs(x[0]);

    for(int i=1; i<n; i += incx){
        double tried = abs(x[i]);
        if(tried < minValue){
            currentMin = i;
            minValue = tried;
        }
    }

    return currentMin+1;

}

CBLAS_INDEX mnblas_icamin(const int n, const void *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    complexe_float_t *X = (complexe_float_t *)x;

    CBLAS_INDEX currentMin = 0;
    float minValue = abs(X[0].real) + abs(X[0].imaginary);

    for(int i=1; i<n; i += incx){
        float tried = abs(X[i].real) + abs(X[i].imaginary);
        if(tried < minValue){
            currentMin = i;
            minValue = tried;
        }
    }

    return currentMin+1;

}

CBLAS_INDEX mnblas_izamin(const int n, const void *x, const int incx){
    if(n<0 || incx <0){
        return 0;
    }

    complexe_double_t *X = (complexe_double_t *)x;

    CBLAS_INDEX currentMin = 0;
    double minValue = abs(X[0].real) + abs(X[0].imaginary);

    for(int i=1; i<n; i += incx){
        double tried = abs(X[i].real) + abs(X[i].imaginary);
        if(tried < minValue){
            currentMin = i;
            minValue = tried;
        }
    }

    return currentMin+1;

}
