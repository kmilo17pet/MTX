/* 
 * File:   mta2d.h
 * Author: juanc_000
 *
 * Created on 5 de octubre de 2014, 01:19 PM
 */

#ifndef MATRIX_MTA_2D_LIBRARY_H
#define	MATRIX_MTA_2D_LIBRARY_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
    
#ifndef  mta_1dnumel   
    #define mta_1dnumel(X)         (sizeof(X)/sizeof(X[0]))                                    /* Number of elements of ARRAY-1D */
#endif
#define MTA_2DPRINT_FORMAT  "% 8.7g\t"
#define mta_2dnrows(X)       mta_1dnumel(X)                                                 /* Number of rows of ARRAY-2D */
#define mta_2dncols(X)       mta_1dnumel(X[0])                                              /* Number of columns of ARRAY-2D */

int mta_2doper(double *out,const double *a,const double *b,const int nf,const int nc,const int nra, const int nca,const int nrb,const int ncb,const int oper);
#define   mta_s_2dadd(OUT,A,B)                  mta_2doper(&OUT[0][0],&A[0][0],&B[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),mta_2dnrows(B),mta_2dncols(B),1)
#define   mta_s_2dsub(OUT,A,B)                  mta_2doper(&OUT[0][0],&A[0][0],&B[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),mta_2dnrows(B),mta_2dncols(B),2)
#define   mta_s_2dpptim(OUT,A,B)                mta_2doper(&OUT[0][0],&A[0][0],&B[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),mta_2dnrows(B),mta_2dncols(B),3)
#define   mta_s_2dppdiv(OUT,A,B)                mta_2doper(&OUT[0][0],&A[0][0],&B[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),mta_2dnrows(B),mta_2dncols(B),4)
#define   mta_s_2dprod(OUT,A,B)                 mta_2doper(&OUT[0][0],&A[0][0],&B[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),mta_2dnrows(B),mta_2dncols(B),5)

int mta_2dsoper(double *out, const double *a,const int nf,const int nc,const int nfa, const int nca, const int oper);
#define   mta_s_2dsadd(OUT,A)                mta_2dsoper(&OUT[0][0],&A[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),1)
#define   mta_s_2dssub(OUT,A)                mta_2dsoper(&OUT[0][0],&A[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),2)
#define   mta_s_2dspptim(OUT,A)              mta_2dsoper(&OUT[0][0],&A[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),3)
#define   mta_s_2dsppdiv(OUT,A)              mta_2dsoper(&OUT[0][0],&A[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),4)
#define   mta_s_2dsprod(OUT,A)               mta_2dsoper(&OUT[0][0],&A[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),5)

int mta_2deye(double *out,const int nf,const int nc);
#define   mta_s_2deye(OUT)                      mta_2deye(&OUT[0][0],mta_2dnrows(OUT),mta_2dncols(OUT))  
int mta_2dpowui(double *out, const double *in,const unsigned int power,const int nfo,const int nco, const int nf, const int nc);
#define   mta_s_2dpowui(OUT_,IN_,POWER_)          mta_2dpowui((&OUT_[0][0]) ,(&IN_[0][0]), POWER_, mta_2dnrows(OUT_), mta_2dncols(OUT_), mta_2dnrows(IN_), mta_2dncols(IN_) )
int mta_2dt(double *out,const double *a,const int nf,const int nc,const int nra, const int nca);
#define  mta_s_2dt(OUT,A)                   mta_2dt(&OUT[0][0],&A[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A))
int mta_2dkoper(double *out,const double *a,const double s,const int nf,const int nc,const int oper);
#define   mta_s_2dkadd(OUT,A,SC)                 mta_2dkoper(&OUT[0][0],&A[0][0],SC,mta_2dnrows(OUT),mta_2dncols(OUT),1)
#define   mta_s_2dksub(OUT,A,SC)                 mta_2dkoper(&OUT[0][0],&A[0][0],SC,mta_2dnrows(OUT),mta_2dncols(OUT),2)
#define   mta_s_2dkprod(OUT,A,SC)                mta_2dkoper(&OUT[0][0],&A[0][0],SC,mta_2dnrows(OUT),mta_2dncols(OUT),3)
#define   mta_s_2dkdiv(OUT,A,SC)                 mta_2dkoper(&OUT[0][0],&A[0][0],SC,mta_2dnrows(OUT),mta_2dncols(OUT),4)
#define   mta_s_2dkpow(OUT,A,SC)                 mta_2dkoper(&OUT[0][0],&A[0][0],SC,mta_2dnrows(OUT),mta_2dncols(OUT),5)
int mta_2dgetsubset(double *outm,const double *m, const int f1, const int f2, const int c1, const int c2,const int nf,const int nc,const int nfm, const int ncm);
#define mta_s_2dgetsubset(OUT,M,F1,F2,C1,C2)   mta_2dgetsubset(&OUT[0][0],&M[0][0],F1,F2,C1,C2,mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M))
#define mta_s_2dgetrow(OUT,M,R)                mta_2dgetsubset(&OUT[0][0],&M[0][0],R,R,0,(mta_2dncols(M)-1),mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M))
#define mta_s_2dgetcol(OUT,M,C)                mta_2dgetsubset(&OUT[0][0],&M[0][0],0,(mta_2dnrows(M)-1),C,C,mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M))
#define mta_s_2dgetrows(OUT,M,F1,F2)           mta_2dgetsubset(&OUT[0][0],&M[0][0],F1,F2,0,(mta_2dncols(M)-1),mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M))
#define mta_s_2dgetcols(OUT,M,C1,C2)           mta_2dgetsubset(&OUT[0][0],&M[0][0],0,(mta_2dnrows(M)-1),C1,C2,mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M))
double mta_2dmax(int *row_m , int *col_m, double *m, const int nf, const int nc);
#define mta_s_2dmax(f_Max,c_Max,MAT)            mta_2dmax(&f_Max,&c_Max,&MAT[0][0],mta_2dnrows(MAT),mta_2dncols(MAT))
double mta_2dmin(int *row_m , int *col_m, double *m, const int nf, const int nc);
#define mta_s_2din(f_Min,c_Min,MAT)            mta_2dmin(&f_Min,&c_Min,&MAT[0][0],mta_2dnrows(MAT),mta_2dncols(MAT))
int mta_2dvcat(double *out, const double *a, const double *b, const int nfo, const int nco, const int nfa, const int nca, const int nfb, const int ncb);
#define mta_s_2dvcat(OUT,A,B)                   mta_2dvcat(&OUT[0][0],&A[0][0],&B[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),mta_2dnrows(B),mta_2dncols(B))
int mta_2dhcat(double *out, const double *a, const double *b, const int nfo, const int nco, const int nfa, const int nca, const int nfb, const int ncb);
#define mta_s_2dhcat(OUT,A,B)                   mta_2dhcat(&OUT[0][0],&A[0][0],&B[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A),mta_2dnrows(B),mta_2dncols(B))
int mta_2dsetsubset(double *out,const double *a, const int f1, const int f2, const int c1, const int c2, const int nfo, const int nco, const int nfa, const int nca);
#define mta_s_2dsetsubset(OUT, A, F1,F2,C1,C2) mta_2dsetsubset(&OUT[0][0],&A[0][0],F1,F2,C1,C2,mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A))
int mtx_2dsetvalsubset(double *out,const double val, const int f1, const int f2, const int c1, const int c2, const int nfo, const int nco);
#define mtx_s_2dsetvalsubset(OUT, VAL, F1,F2,C1,C2) mtx_2dsetvalsubset(&OUT[0][0],(VAL),F1,F2,C1,C2,mta_2dnrows(OUT),mta_2dncols(OUT))
int mta_2dcpy(double *out,const double *a,const int nf,const int nc,const int nfa, const int nca);
#define   mta_s_2dcpy(OUT,A)                    mta_2dcpy(&OUT[0][0],&A[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(A),mta_2dncols(A)) 
int mta_2dvecgen(double *out,const double init,const double inc, const double endi, const int row, const int col, const int nf, const int nc);
#define mta_s_2dvecgen(OUT,INIT,INCR,END,ROW,COL) mta_2dvecgen(&OUT[0][0],INIT,INCR,END,ROW,COL,mta_2dnrows(OUT),mta_2dncols(OUT))        
void mta_2ddisp(const double *array2d,int nf,int nc);      
#define mta_s_2ddisp(M)   mta_2ddisp(&M[0][0],mta_2dnrows(M),mta_2dncols(M))   
#define mta_2dshow(M)       printf(""#M"[%d][%d]=",mta_2dnrows(M),mta_2dncols(M)); \
                            mta_s_2ddisp(M)                                      

int mta_2dgetsubsetbyidx(double *out, const double *m, const int *mf, const int *mc , const int nf, const int nc, const int nfm, const int ncm, int lmf, int lmc);
#define mta_s_2dgetsubsetbyidx(OUT,M,FM,CM)      mta_2dgetsubsetbyidx(&OUT[0][0],&M[0][0],&FM[0],&CM[0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M),mta_1dnumel(FM),mta_1dnumel(CM))
#define mta_2dgetcolsbyidx(OUT,M,CM)           mta_2dgetsubsetbyidx(&OUT[0][0],&M[0][0],&((int *)0)[0],&CM[0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M),-1,mta_1dnumel(CM))
#define mta_2dgetrowsbyidx(OUT,M,FM)           mta_2dgetsubsetbyidx(&OUT[0][0],&M[0][0],&FM[0],&((int *)0)[0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M),mta_1dnumel(FM),-1)

int mta_2dinv(double *out, const double *m, double *perm, const int nfout, const int ncout, const int nfm, const int ncm, int nfperm, int ncperm);
#define mta_s_2dinv(OUT,M,PERM)             mta_2dinv(&OUT[0][0],&M[0][0],&PERM[0][0],mta_2dnrows(OUT),mta_2dncols(OUT),mta_2dnrows(M),mta_2dncols(M),mta_2dnrows(PERM),mta_2dncols(PERM))   

void mta_2dinit(double *ar,const int nf,const int nc,const double val);
#define mta_s_2dinit(A,val)    mta_2dinit(&A[0][0],mta_2dnrows(A),mta_2dncols(A),val)
int mta_2dswaprows(double *ar, const int r1, const int r2, const int nf, const int nc);
#define mta_s_2dswaprows(A,R1,R2)       mta_2dswaprows(&A[0][0], R1, R2, mta_2dnrows(A), mta_2dncols(A))
int mta_2dswapcols(double *ar, const int c1, const int c2, const int nf, const int nc);
#define mta_s_2dswapcols(A,C1,C2)       mta_2dswapcols(&A[0][0], C1, C2, mta_2dnrows(A), mta_2dncols(A))

#ifdef	__cplusplus
}
#endif

#endif	/* MTA2D_H */

