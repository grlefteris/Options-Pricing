#pragma once
#include <cmath>
#include "matrix.h"

matrix cholesky(const matrix &A) {
    const size_t n = A.rows();
    matrix L(n,n);

    L(0,0) = sqrt( A(0,0) );

    for (size_t i = 1; i < n; i++) {
        L(i,0) = A(i,0) / L(0,0);
    }

    double sum;
    for (size_t i = 1; i < n-1; i++) {
        sum = 0.0;
        for (size_t k = 0; k < i; k++) {
            sum += L(i,k) * L(i,k);
        }
        L(i,i) = sqrt( A(i,i) - sum );
        for (size_t j = i+1; j < n; j++) {
            sum = 0.0;
            for (size_t k = 0; k < i; k++) {
                sum += L(j,k) * L(i,k);
            }
            L(j,i) = ( A(j,i) - sum ) / L(i,i);
        }
    }
    sum=0;

    for (size_t k = 0; k < n-1; k++) {
        sum += L(n-1,k) * L(n-1,k);
    }
    L(n-1,n-1) = sqrt( A(n-1,n-1) - sum );

    return L;
}
