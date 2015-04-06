/* 
 * File:   mta.h
 * Author: juanc_000
 *
 * Created on 5 de octubre de 2014, 12:54 PM
 */

#ifndef MATRIX_MTA_1D_LIBRARY_H
#define	MATRIX_MTA_1D_LIBRARY_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
    
#define mta_1dnumel(X)         (sizeof(X)/sizeof(X[0]))                                    /* Number of elements of ARRAY-1D */

/*Utiliy functions for 1D-arrays-------------------------------------------------------------------------------------------------*/
#define mta_s_reg_update(X,Y)  mta_update_reg(X,Y,mta_1dnumel(X))                        /* Regression vector update */
void mta_update_reg(double regarray[],const double newp,const int length);
#define mta_s_medfilt(X,Y)    mta_median_filt(X,Y,mta_1dnumel(X))                        /* Median Filter Point by Point */
double mta_median_filt(double regarray[],const double newp,const int length);

void mta_1dinit(double ar[],int arl,double val);
#define mta_s_1dinit(A,v)         mta_1dinit(A,mta_1dnumel(A),v)
double mta_1dpolyval(double ar[],const int arl,const double val);
double mta_1dsum(const double ar[],const int arl);
double mta_1dprod(const double ar[],const int arl);
#define mta_s_1dsum(A)                  mta_1dsum(A,mta_1dnumel(A))
#define mta_s_1dmean(A)                 (mta_s_1dsum(A)/mta_1dnumel(A))
#define mta_s_1dpolyval(A,val)          mta_1dpolyval(A,mta_1dnumel(A),val)
#define mta_s_1dprod(A)                 mta_1dprod(A,mta_1dnumel(A))

int mta_1doper(double out[],const double a[],const double b[],const int lout,const int la,const int lb,const int oper);

#define MTA_ADD     1
#define MTA_SUB     2
#define MTA_TIM     3
#define MTA_DIV     4

#define mta_s_1dadd(OUT,A,B)    mta_1doper(OUT,A,B,mta_1dnumel(OUT),mta_1dnumel(A),mta_1dnumel(B),MTA_ADD)
#define mta_s_1dsub(OUT,A,B)    mta_1doper(OUT,A,B,mta_1dnumel(OUT),mta_1dnumel(A),mta_1dnumel(B),MTA_SUB)
#define mta_s_1dtim(OUT,A,B)    mta_1doper(OUT,A,B,mta_1dnumel(OUT),mta_1dnumel(A),mta_1dnumel(B),MTA_TIM)
#define mta_s_1ddiv(OUT,A,B)    mta_1doper(OUT,A,B,mta_1dnumel(OUT),mta_1dnumel(A),mta_1dnumel(B),MTA_DIV)
void mta_1dminus_one(double ar[],const int l);
#define mta_s_1dn(A)      mta_1dminus_one(A,numel_1d(A))

#define mta_s_1dfx(OUT, A, fx)  mta_1dfx(OUT, A, mta_1dnumel(OUT), mta_1dnumel(A), fx)
void mta_1dfx(double out[], double ar[], const int lout, const int la, double (*fx)(double));

void mta_1dvecgen(double *out,const double init,const double inc, const double endi, const int n);
#define array1d_vec_gen(OUT,INIT,INCR,END)      mta_1dvecgen(&OUT[0], INIT, INCR , END, mta_1dnumel(OUT))
 
int mta_1dmin(double ar[], int la);
int mta_1dmax(double ar[], int la);
#define mta_s_1dmin(A)      mta_1dmin(A, mta_1dnumel(A)) 
#define mta_s_1dmax(A)      mta_1dmax(A, mta_1dnumel(A)) 

#ifdef	__cplusplus
}
#endif

#endif	/* MTA_H */

