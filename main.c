/* 
 * Evaluate the expresion X= (A*B)/(alpha + B'*A*B) 
 *  A is a 4-square random matrix and B is a 4-column random matrix
 * alpha is scalar
 */
#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"

int main(void){
    /*Define the  using the matrix notation*/
    matrix A = mtx_rand(4,4);
    matrix B = mtx_rand(4,1);
    matrix X = NULL; //holds the expression
    matrix AB = NULL; //holds the auxiliar A*B result
    matrix BtAB = mtx_new(1,1); // holds the auxiliar B'*A*B result 
    double alpha = 0.8;
    
    /*Compute the operations*/
    AB = mtx_prod(1.0, A,B); // AB = A*B
    mtx_dgemm(1.0, B, 't', AB, '.', 0, BtAB); // compute B'*A*B using the generalized matrix product routine mtx_dgemm
    X = mtx_koper(AB, '/', (alpha + BtAB->pos[0][0]) );
    /*Display the matrices*/
    mtx_dispn(A,B,AB,BtAB,X);
    /*Release the heap*/
    mtx_del(A);
    mtx_del(B);
    mtx_del(AB);
    mtx_del(BtAB);
    mtx_del(X);
    return EXIT_SUCCESS;
}