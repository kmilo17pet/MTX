/*
 * MTX, a BLAS library (level 3) for ANSI-C  v3.4
 *  	Copyright (C) 2013 Eng. Juan Camilo Gómez C. MSc. (kmilo17pet@gmail.com)
 *
 *	MTX is free software: you can redistribute it and/or modify it
 *	under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	MTX is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef MATRIX_MTX_LIBRARY_H
#define	MATRIX_MTX_LIBRARY_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <stdarg.h>
#include <stdint.h>

    
struct {
    unsigned char rows,cols,mem;
    double **pos; 
}
#if defined(__GNUC__)
typedef  _Matrix, *matrix;
#elif defined(_MSC_VER)
typedef struct _Matrix *matrix;
#else
typedef _Matrix, *matrix;
#endif
    
    
#define MTX_PRINT_FORMAT    "%8.7g\t"
#define MATRIX              matrix
#define Matrix              matrix
#define mtx_show(M)         fputs("["#M"]",stdout);mtx_disp(M);\
         
#define mtx_dispn(args...)  mtx_ndisp(50,args,-1)
matrix mtx_ndisp(int n, ...);

#define mtxdef(_VAR_)       matrix _VAR_=NULL
#define mtx_del(M)          _mtx_del(M);M=NULL
#define mtx_length(M)       ( (M==NULL)? -1 : (((M->rows)>(M->cols))? M->rows : M->cols) )      /*Length of Matrix X - max(rows,cols)*/    
#define mtx_isempty(M)      ( (M==NULL)? 1 : ((M->cols<=0)|(M->rows<=0)) )                      /* TRUE if M=[] */
#define mtx_isrow(M)        ( (M==NULL)? 0 : (((M->cols)>0)&&((M->rows)==1)) )                  /* TRUE if size(M)=1xN */
#define mtx_iscolumn(M)     ( (M==NULL)? 0 : (((M->rows)>0)&&((M->cols)==1)) )                  /* TRUE if size(M)=Nx1 */
#define mtx_isvector(M)     ( (M==NULL)? 0 : (mtx_isrow(M) || mtx_iscolumn(M)) )               /* TRUE if M is vector (row or column)*/
#define mtx_numel(M)        ( (M==NULL)? -1 :  (((M)->rows)*((M)->cols)) )
    
#define mtx_issquare(M)     ( (M==NULL)? 0 :  ((M->rows)==(M->cols)) )
#define mtx_ismultiply(A,B) ( (M==NULL)? 0 : ((A->cols) == (B->rows)) )
#define mtx_samedim(A,B)    ( (M==NULL)? 0 : (((A->rows) == (B->rows)) && ((A->cols) != (B->cols))) )
#define mtx_lastc(M)        ((M==NULL)? -1 :(M->cols)-1)
#define mtx_lastr(M)        ((M==NULL)? -1 :(M->rows)-1)

#define mtx_setrow(M,A,R)              mtx_setsubset(M,A,R,R,0,mtx_lastc(M))                      /* M(R,:)=A */
#define mtx_setcol(M,A,C)              mtx_setsubset(M,A,0,mtx_lastr(M),C,C)                      /* M(R,:)=A */

#define mtx_getrow(M,R)                mtx_getsubset(M,R,R,0,mtx_lastc(M))                        /* M(R,:) Get column from matrix */
#define mtx_getcol(M,C)                mtx_getsubset(M,0,mtx_lastr(M),C,C)                        /* M(:,C) Get row from matrix */
#define mtx_getrows(M,F1,F2)           mtx_getsubset(M,F1,F2,0,mtx_lastc(M))                      /* M(F1:F2,:) Get columns from matrix */
#define mtx_getcols(M,C1,C2)           mtx_getsubset(M,0,mtx_lastr(M),C1,C2)                      /* M(:,C1:C2) Get rows from matrix */    



matrix mtx_new(const unsigned char rows,const unsigned char cols);
void _mtx_del(const matrix M);
matrix mtx_cpy(const matrix M);
matrix mtx_eye(unsigned char n, const double alpha);
matrix mtx_diag(const matrix m);
double mtx_trace(const matrix A);
void mtx_disp(const matrix M);
matrix mtx_t(const matrix A);                   
matrix mtx_gadd(const double alpha, const matrix A, const double beta, const matrix B);
matrix mtx_ptpprod(const matrix A, const matrix B);
matrix mtx_ptpdiv(const matrix A, const matrix B);
matrix mtx_koper(const matrix A, const char oper, const double k);
matrix mtx_rand(unsigned rows, unsigned char cols);
matrix mtx_prod(const double alpha, const matrix A, const matrix B);

int mtx_OUT_equal_AxB(matrix OUT, const double alpha, const matrix A, const matrix B);
int mtx_A_equal_A_plus_B(matrix A, const matrix B);
int mtx_A_equal_A_sub_B(matrix A, const matrix B);
int mtx_A_equal_kA(matrix A, const double k);
int mtx_OUT_equal_kA(matrix OUT,const matrix A, const double k);
matrix mtx_inv(const matrix X);
matrix mtx_linsolve(const matrix A, const matrix B);
matrix mtx_hcat(const matrix m1, const matrix m2);
matrix mtx_vcat(const matrix m1, const matrix m2);
double mtx_lu(matrix L, matrix U, const matrix M);
double mtx_det(const matrix M);
int mtx_swaprows(matrix A, int r1, int r2);
int mtx_swapcols(matrix A, int c1, int c2);
matrix mtx_powui(const matrix m, unsigned int power); //not working yet
int mtx_memcpy(matrix dst, matrix scr);

matrix mtx_mean(const matrix m);
double mtx_cumsum(const matrix m);
matrix mtx_rpinv(const matrix A, double tol);
matrix mtx_lpinv(const matrix A, double tol);
double mtx_cumprod(const matrix m);
matrix mtx_produ(const matrix m);
matrix mtx_sum(const matrix m);
double mtx_cummax(const matrix m);
matrix mtx_max(const matrix m);
double mtx_cummin(const matrix m);
matrix mtx_min(const matrix m);

matrix mtx_2d2mtx(const double *array2d,const int nf,const int nc);
#define mtx_2dtomtx(A)      mtx_2d2mtx(&A[0][0],(sizeof(A)/sizeof(A[0])), (sizeof(A[0])/sizeof(A[0][0])))
matrix mtx_1d2mtx(const double *array1d, const int arraylength);
#define mtx_1dtomtx(A)      mtx_1d2mtx(A,(sizeof(A)/sizeof(A[0])))

double mtx_colspprod(const matrix A, const int j, const int k);
matrix mtx_grams(const matrix M, matrix R); // -> ok
#define mtx_qr(A,R)     mtx_grams(A,R)

double mtx_norm_inf(const matrix M);
matrix mtx_expm(const matrix M, const double alpha); // exponential matrix (using pade aprox) -> ok
int mtx_dgemm(const double alpha, const matrix a, char transa, const matrix b, char transb, const double beta, matrix c); //computes C = alpha * A * B and related operations.
int mtx_dgema(const double alpha, const matrix a, char transa, const double beta, const matrix b, char transb); // under construction
matrix mtx_fxptp(const matrix A, double (*fx)(double));

#define mtx_vcatn(args...)  mtx_nvcat(50,args,NULL)
matrix mtx_nvcat(int n, ...);
#define mtx_hcatn(args...)  mtx_nhcat(50,args,NULL)
matrix mtx_nhcat(int n, ...);

matrix mtx_lspcf(double *X, double *Y, int n, int m);
matrix mtx_dot(matrix A, matrix B);
matrix mtx_cholesky(matrix A);
        
matrix mtx_getsubset(const matrix m, int f1, int f2, int c1, int c2);
void mtx_setsubset(matrix m1, const matrix m2,int f1,int f2,int c1,int c2);


#ifdef	__cplusplus
}
#endif

#endif	/* MATRIX_MTX_LIBRARY_H */

