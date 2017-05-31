#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"

int main(void){
    mtx_assign_heap_wrappers(calloc, free);
    double a[][4]= {
     1,     0,     0,     0,
     0,     0,     0,     1,
     5,     0,     0,     1,
     1,     2,     3,     4,
                    };
    matrix A = mtx_2dtomtx(a);
    matrix U = mtx_new(A->rows, A->cols);
    matrix S = mtx_new(A->cols, A->cols);
    matrix V = mtx_new(A->cols, A->cols);
    mtx_svd(A,U,S,V);
    matrix t1 = mtx_prod(1,mtx_prod(1,U,S),mtx_t(V));
    #ifdef MTX_PRINTOUT
    mtx_dispn(A,U,S,V,t1);
    #endif
    mtx_del(A);
    mtx_del(U);
    mtx_del(S);
    mtx_del(V);
    mtx_del(t1);
    return EXIT_SUCCESS;
}