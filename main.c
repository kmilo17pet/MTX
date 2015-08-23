#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"

int main(void){    
    double a[3][3]={
    {0,1,8},
    {3,0,2},
    {0,1,1},
    };
    double b[3][3]={
    {10,11,12},
    {13,14,15},
    {16,17,19},
    };
    double c[3][3]={
    {1,2,3},
    {4,5,6},
    {7,8,9},
    };
    double z[5][5]={
        {1,    1,    1,    1,    1},
        {1,    2,    3,    4,    5},
        {1,    3,    6,   10,   15},
        {1,    4,   10,   20,   35},
        {1,    5,   15,   35,   70},
    };
    double y[6]={8, 8 , 9, 1 , 2 ,3};
    matrix Y= mtx_1dtomtx(y);
    matrix W = mtx_2dtomtx(z);
    z[0][0]=.2;
    y[2]=-2;
    matrix A=mtx_2dtomtx(a);

    matrix Z=mtx_2dtomtx(z);
    matrix R = mtx_new(A->rows, A->cols);
    matrix Q = mtx_qr(A, R);
    mtx_dispn(A,Q,R,mtx_prod(1.0,Q,R),mtx_inv(A),mtx_inv(Z), W, Y);  
    long i;
    matrix L=NULL;
    matrix U=mtx_rand(100,100);
    while(1){
        L=mtx_inv(U);
        mtx_del(L);
        puts("running");
    }
    mtx_del(Y);
    mtx_del(W);
    mtx_del(R);
    mtx_del(A);
    mtx_del(Q);mtx_del(Z);mtx_del(U);
    
    return EXIT_SUCCESS;
}