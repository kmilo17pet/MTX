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


double x[]={    0,     0.1000,    0.2000,    0.3000,    0.4000,    0.5000,    0.6000,    0.7000,    0.8000,    0.9000,    1.0000,    1.1000,    1.2000,    1.3000,
            1.4000,    1.5000,    1.6000,    1.7000,    1.8000,    1.9000,    2.0000,    2.1000,    2.2000,    2.3000,    2.4000,    2.5000,    2.6000,    2.7000,
            2.8000,    2.9000,    3.0000,    3.1000,    3.2000,    3.3000,    3.4000,    3.5000,    3.6000,    3.7000,    3.8000,    3.9000,    4.0000,    4.1000,
            4.2000,    4.3000,    4.4000,    4.5000,    4.6000,    4.7000,    4.8000,    4.9000,    5.0000,};
double y[]={1.2000,    1.2530,    1.3120,    1.3770,    1.4480,    1.5250,    1.6080,    1.6970,    1.7920,    1.8930,    2.0000,    2.1130,    2.2320,    2.3570,
            2.4880,    2.6250,    2.7680,    2.9170,    3.0720,    3.2330,    3.4000,    3.5730,    3.7520,    3.9370,    4.1280,    4.3250,    4.5280,    4.7370,
            4.9520,    5.1730,    5.4000,    5.6330,    5.8720,    6.1170,    6.3680,    6.6250,    6.8880,    7.1570,    7.4320,    7.7130,    8.0000,    8.2930,
            8.5920,    8.8970,    9.2080,    9.5250,    9.8480,   10.1770,   10.5120,   10.8530,   11.2000,};

int main(void){
    double a[][4]= {
                    1,0,0,0,
                    0,1,0,0,
                    0,0,1,0,
                    0,0,2,2,
                    };
    matrix A = mtx_2dtomtx(a);
    matrix U = mtx_new(4,4);
    matrix S = mtx_new(4,4);
    matrix V = mtx_new(4,4);
    printf("\r\n svd retval = %d\r\n",mtx_svd(A,U,S,V));
    matrix t1 = mtx_prod(1,mtx_prod(1,U,S),mtx_t(V));
    mtx_dispn(A,U,S,V,t1);
    return EXIT_SUCCESS;
}