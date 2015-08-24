/* solving the matrix equation A*x=b using MTX */
#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"

int main(void){
    double a[3][3]= {
                    3.1, 1.3,-5.7,
                    1.0, -6.9, 5.8,
                    3.4, 7.2, -8.8
                    };
    double b[2][2]= {
                    0.0357,    0.9340,
                    0.8491,    0.6787 
                    };
    double c[3][2]={
                    0.7577,    0.6555,
                    0.7431,    0.1712,
                    0.3922,    0.7060,
                    };
    double d[3][1]={
                    1,
                    2,
                    3
                    };
    matrix A = mtx_2dtomtx(a);
    matrix B = mtx_2dtomtx(b);
    matrix C = mtx_2dtomtx(c);
    matrix D = mtx_2dtomtx(d);
    matrix X = NULL;
    matrix Y = NULL;
    do{
        X = mtx_sylvester(A,B,C);
        mtx_del(X);
    }while(1);
    mtx_del(A);
    mtx_del(B);
    mtx_del(C);
    return EXIT_SUCCESS;
}