#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"

int main(void){    
    double a[3][3]={
    {1,2,0},
    {0,-5,6},
    {7,0,9},
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
    matrix A=mtx_2dtomtx(a);
    matrix B=mtx_2dtomtx(b);
    matrix C=mtx_2dtomtx(c);
    matrix Z=NULL;
        mtx_dispn(A,B,C);

    mtx_dgemm(-0.01, A,'-', B,'-', 1.0,C);
        mtx_dispn(A,B,C, mtx_expm(A, 0.01));
        


    //matrix B=mtx_2dtomtx(b);
    //matrix C=mtx_2dtomtx(c);
    //matrix Z=mtx_2dtomtx(z);
    //Z=mtx_inv(A);
    //mtx_dispn(A,Z);
    //mtx_dispn(A,B,C,Z);
    //mtx_dispn(Z,mtx_cholesky(Z));
    
    return EXIT_SUCCESS;
}