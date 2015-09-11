/* 
 * Evaluate the expresion X= (A*B)/(alpha + B'*A*B) 
 *  A is a 4-square random matrix and B is a 4-column random matrix
 * alpha is scalar
 */
#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"
#include <time.h>
#include <unistd.h>
#include <sys/time.h>

struct timeval  tv1, tv2;
clock_t begin, end;
double time_spent;

int main(void){
    double v[][5]= {
     0,     0, 5   ,1, 3,
     0,    5,   -5  ,1, 3,
     1,     0,   5,  1, 0,
                    };
    double a[][3]={
    0,0,0,
    0,9,0,
    3,0,0,
    };
    matrix A = mtx_2dtomtx(v);
    matrix U = mtx_new(A->rows, A->cols);
    matrix S = mtx_new(A->cols, A->cols);
    matrix V = mtx_new(A->cols, A->cols);
    mtx_svd(A,U,S,V);
    matrix t1 = mtx_prod(1,mtx_prod(1,U,S),mtx_t(V));
    mtx_dispn(A,U,S,V,t1);
    return EXIT_SUCCESS;
}