#include "mta2d.h"

/*============================================================================*/
int mta_2doper(double *out,const double *a,const double *b,const int nf,const int nc,const int nra, const int nca,const int nrb,const int ncb,const int oper){
   int i,j,y;
   switch (oper){
   case 1:
       if ( ((nra!=nrb) || ((nf!=nra)||(nf!=nrb)))   ||   ((nca!=ncb) || ((nc!=nca)||(nc!=ncb)))  ) return -1;
       for (i=0;i<nf;i++){
                 for (j=0;j<nc;j++){
                 out[i*nc+j]=a[i*nc+j]+b[i*nc+j];
        }}
       break;
   case 2:
       if ( ((nra!=nrb) || ((nf!=nra)||(nf!=nrb)))   ||   ((nca!=ncb) || ((nc!=nca)||(nc!=ncb)))  ) return -1;
       for (i=0;i<nf;i++){
                 for (j=0;j<nc;j++){
                 out[i*nc+j]=a[i*nc+j]-b[i*nc+j];
        }}
       break;
   case 3:
       if ( ((nra!=nrb) || ((nf!=nra)||(nf!=nrb)))   ||   ((nca!=ncb) || ((nc!=nca)||(nc!=ncb)))  ) return -1;
       for (i=0;i<nf;i++){
                 for (j=0;j<nc;j++){
                 out[i*nc+j]=a[i*nc+j]*b[i*nc+j];
       }}
       break;
   case 4:
       if ( ((nra!=nrb) || ((nf!=nra)||(nf!=nrb)))   ||   ((nca!=ncb) || ((nc!=nca)||(nc!=ncb)))  ) return -1; 
       for (i=0;i<nf;i++){
                 for (j=0;j<nc;j++){
                 out[i*nc+j]=a[i*nc+j]/b[i*nc+j];
       }}
       break;
   case 5:
       if(  ((nca != nrb)||(nc!=ncb)) || (nf!=nra)  ) return -1;
        for(i =0;i<nra; i++){
            for(j=0;j<ncb; j++){
                out[i*nc+j]=0.0;
                for(y=0;y<nca; y++){
                    out[i*nc+j]=out[i*nc+j]+a[i*nca+y]*b[y*ncb+j]; 
                }
            }
        }
       break;
   default:
       if ( ((nra!=nrb) || ((nf!=nra)||(nf!=nrb)))   ||   ((nca!=ncb) || ((nc!=nca)||(nc!=ncb)))  ) return -1;
       for (i=0;i<nf;i++){
                 for (j=0;j<nc;j++){
                 out[i*nc+j]=a[i*nc+j]+b[i*nc+j];
        }}
       break;
   }
   return 0;
}
/*============================================================================*/
int mta_2dsoper(double *out, const double *a,const int nf,const int nc,const int nfa, const int nca, const int oper){ //+-*/x
    if ( nf!=nfa && nc!=nca  ) return -1;
    int i,j,y;
    int nnf=nf,nnc=nc;
    if (oper!=5) {nnf=0;nnc=0;}
    
    double temp[nnf][nnc];
    
    switch (oper){
        case 1:
            for (i=0;i<nf;i++){
                for (j=0;j<nc;j++){
                out[i*nc+j]+=a[i*nc+j];
            }}
            break;
        case 2:   
            for (i=0;i<nf;i++){
                for (j=0;j<nc;j++){
                out[i*nc+j]-=a[i*nc+j];
            }}
            break;
        case 3:
            for (i=0;i<nf;i++){
                for (j=0;j<nc;j++){
                out[i*nc+j]*=a[i*nc+j];
            }} 
            break;
        case 4:
            for (i=0;i<nf;i++){
                for (j=0;j<nc;j++){
                out[i*nc+j]/=a[i*nc+j];
            }}    
            break;
        case 5:
            for(i=0 ;i<nf; i++){
                for(j=0; j<nc; j++){
                    temp[i][j]=0.0;
                    for(y=0;y<nc; y++){
                        temp[i][j]=temp[i][j]+out[i*nc+y]*a[y*nc+j]; //_m(c,i,j)=_m(c,i,j)+_m(m1,i,y)*_m(m2,y,j);
                    }
                }
            }
           for(i=0 ;i<nf; i++){
                for(j=0; j<nc; j++){
                    out[i*nc+j]=temp[i][j];
           }}
           break;            
        default:
            break;
    }
    return 0;
}
/*============================================================================*/
int mta_2dpowui(double *out, const double *in,const unsigned int power,const int nfo,const int nco, const int nf, const int nc){ //not working yet
    if( nf!=nc || nfo!=nco || nf!=nfo) return -1;
    mta_2deye(out,nf,nc);
    if (power==0)  return;    
    int k;
    for (k=0;k<power;k++)
        mta_2dsoper(out, in,nf,nf,nf, nf,  5);
    return 0;
}
/*============================================================================*/
int mta_2deye(double *out,const int nf,const int nc){
    if (nf<=0 || nc<=0 || nf!=nc) return -1;
    int i,j;
    for (i=0;i<nf;i++){
                 for (j=0;j<nc;j++){
                 out[i*nc+j]= (i==j)? 1.0 : 0.0;
    }}
    return 0;
}
/*============================================================================*/
int mta_2dt(double *out,const double *a,const int nf,const int nc,const int nra, const int nca){
    if ( nf!=nca || nc!=nra ) return -1;
    int i,j;
    for (i=0;i<nra;i++){
        for (j=0;j<nca;j++){
        out[j*nc+i]=a[i*nca+j];
    }}
    return 0;
}
/*============================================================================*/
int mta_2dkoper(double *out,const double *a,const double s,const int nf,const int nc,const int oper){
   int i,j;
   switch (oper){
    case 1:
       for (i=0;i<nf;i++){
                 for (j=0;j<nc;j++){
                 out[i*nc+j]=a[i*nc+j]+s; 
        }}
       break;
    case 2:
        for (i=0;i<nf;i++){
                  for (j=0;j<nc;j++){
                  out[i*nc+j]=a[i*nc+j]-s;
         }}
        break;
    case 3:
        for (i=0;i<nf;i++){
                  for (j=0;j<nc;j++){
                  out[i*nc+j]=a[i*nc+j]*s;
        }}
        break;
    case 4:
        for (i=0;i<nf;i++){
                  for (j=0;j<nc;j++){
                  out[i*nc+j]=a[i*nc+j]/s;
        }}
        break;
    case 5:
        for (i=0;i<nf;i++){
                  for (j=0;j<nc;j++){
                  out[i*nc+j]=pow(a[i*nc+j],s);
        }}
        break;
    default:
        for (i=0;i<nf;i++){
                  for (j=0;j<nc;j++){
                  out[i*nc+j]=a[i*nc+j]+s;
        }}
        break;
    }
   return 0;
}
/*============================================================================*/
int mta_2dgetsubset(double *out,const double *m, const int f1, const int f2, const int c1, const int c2,const int nf,const int nc,const int nfm, const int ncm){
    int ft=f2-f1+1;
    int ct=c2-c1+1;
    int f,c;
    if (ft<=0) ft=1;
    if (ct<=0) ct=1;
    for (f=f1;f<=f2;f++){
        for (c=c1;c<=c2;c++)
            out[(f-f1)*nc+(c-c1)]=m[f*ncm+c];
    }  
    return 0;
}
/*============================================================================*/
double mta_2dmax(int *row_m , int *col_m, double *m, const int nf, const int nc){
    int f,c;
    double maxval,temp;
    (*row_m)=0;
    (*col_m)=0;
    maxval=m[0];
    for (f=0;f<nf;f++){
        for (c=0;c<nc;c++){
            temp=m[f*nc+c];
            if (temp>maxval){
                (*row_m)=f;
                (*col_m)=c;
                maxval=temp;
            }
        }   
    }
    return maxval;
}
/*============================================================================*/
double mta_2dmin(int *row_m , int *col_m, double *m, const int nf, const int nc){
    int f,c;
    double minval,temp;
    (*row_m)=0;
    (*col_m)=0;
    minval=m[0];
    for (f=0;f<nf;f++){
        for (c=0;c<nc;c++){
            temp=m[f*nc+c];
            if (temp<minval){
                (*row_m)=f;
                (*col_m)=c;
                minval=temp;
            }
        }   
    }
    return minval;
}
/*============================================================================*/
int mta_2dhcat(double *out, const double *a, const double *b, const int nfo, const int nco, const int nfa, const int nca, const int nfb, const int ncb){
    if ( nfo!=nfa || nfb!=nfa || nfb!=nfo ) return -1;
    int f,c,c1=nca;
    int nc=nca+ncb;
    for (f=0;f<nfa;f++){
        for (c=0;c<nc;c++){
            out[f*nco + c] = (c>=c1)? b[f*ncb+(c-c1)] : a[f*nca+c];
        }
    }
    return 0;
}
/*============================================================================*/
int mta_2dvcat(double *out, const double *a, const double *b, const int nfo, const int nco, const int nfa, const int nca, const int nfb, const int ncb){
    if ( nco!=nca || ncb!=nca || ncb!=nco ) return -1;
    int f,c,f1=nfa;
    int nf=nfa+nfb;
    for (f=0;f<nf;f++){
        for (c=0;c<nca;c++){
            out[f*nco + c] = (f>=f1)? b[(f-f1)*ncb+c] : a[f*nca+c];
        }
    }
    return 0;
}
/*============================================================================*/
int mta_2dsetsubset(double *out,const double *a, const int f1, const int f2, const int c1, const int c2, const int nfo, const int nco, const int nfa, const int nca){
    if ( (f2-f1)!=(nfa-1) || (c2-c1)!=(nca-1) || nfa>nfo || nca>nco) return -1;
    int f,c;
    for (f=f1;f<=f2;f++){
        for (c=c1;c<=c2;c++){
            out[f*nco+c]=a[(f-f1)*nca + (c-c1)];
        }
    }
    return 0;
}
/*============================================================================*/
int mtx_2dsetvalsubset(double *out,const double val, const int f1, const int f2, const int c1, const int c2, const int nfo, const int nco){
    int f,c;
    for (f=f1;f<=f2;f++){
        for (c=c1;c<=c2;c++){
            out[f*nco+c]=val;
        }
    }
    return 0;
}
/*============================================================================*/
int mta_2dcpy(double *out,const double *a,const int nf,const int nc,const int nfa, const int nca){
    if ( nf!=nfa && nc!=nca ) return -1;
    int i,j;
    for (i=0;i<nfa;i++){
        for (j=0;j<nca;j++){
            out[i*nc+j]=a[i*nca+j];
    }}
    return 0;
}
/*===========================================================================================================*/
int mta_2dvecgen(double *out,const double init,const double inc, const double endi, const int row, const int col, const int nf, const int nc){
    int k;
    double total=((endi-init)/inc)+1;
    if( nf<total || nc<total) return -1;
    double val=init;   
    if (total<0) total=1;
    for (k=0;k<total;k++){
        if (row>=0) out[row*nc+k]=val;
        if (col>=0) out[k*nc+col]=val;
        val+=inc;
    }
    return 0;
}
/*============================================================================*/
void mta_2ddisp(const double *array2d,int nf,int nc){
    double mvalue;
    int i,j;
    putchar('\n');
    for (i=0;i<nf;i++){
        for (j=0;j<nc;j++){
            mvalue=array2d[i*nc+j];
            if (fabs(mvalue)<1e-7) mvalue=0;
            printf(MTA_2DPRINT_FORMAT,mvalue);
        }
    putchar('\n');
    }
    putchar('\n');
}
/*============================================================================*/
int mta_getsubsetbyidx(double *out, const double *m, const int *mf, const int *mc , const int nf, const int nc, const int nfm, const int ncm, int lmf, int lmc){
    int f,c,i,j,itf,itc;
    int ii,jj;

    itf=(lmf<0)? (nf): lmf;
    itc=(lmc<0)? (nc): lmc;

    i=0;j=0;
    ii=0;jj=0;
    while(ii<itf){
          f=(lmf<0)? i : mf[ii];
          ii++;
          if (f<0) continue;
          j=0;jj=0;
          while(jj<itc){
              c=(lmc<0)? j : mc[jj];
              jj++;
              if (c<0) continue;
              out[i*nc+j]=m[f*ncm+c];
              j++;
          }
          i++;
    }
    return 0;
}
/*===========================================================================================================*/
void mta_2dinit(double *ar,const int nf,const int nc,const double val){
    int i,j;
    for (i=0;i<nf;i++){
        for (j=0;j<nc;j++){
            ar[i*nc+j]=val;
    }}
}
/*============================================================================*/       
int mta_2dinv(double *C, const double *X, double *A, const int nfc, const int ncc, const int nfx, const int ncx, int nfa, int nca){
    if(nfc!=ncc || nfx!=ncx || nfa!=nca || nfc!=nfx || nfc!=nfa || nfa!=nfx) return -1;
    mta_2deye(C, nfc, ncc);
    mta_2dcpy(A, X, nfa, nca, nfx, ncx); 
    int k,j; 
    int l,p;   
    double y,m,n;
    for(j=0; j<nca; j++){
        p=j;
        while(fabs(A[p*nca+j]) < 1E-15){ 
            p++;
            if (p >= nfa){
                return -1; 
            }
        }
        y=A[p*nca+j]; 
        for(k=0; k<nca; k++){
            m=A[p*nca+k]; 
            A[p*nca+k]=A[j*nca+k];
            A[j*nca+k]=m/y; 
            n=C[p*nca+k]; 
            C[p*nca+k]=C[j*nca+k];   
            C[j*nca+k]=n/y;
        }
        
        for(l=j+1; l<nca; l++){
            m=A[l*nca+j]; 
            for(k=0; k<nca; k++){
                A[l*nca+k] = A[l*nca+k] - m*A[j*nca+k];
                C[l*nca+k] = C[l*nca+k] - m*C[j*nca+k];
            }
        }
    }
    
    for (j=nca-1; j>=0; j--) 
        for (k=j-1; k>=0; k--){
            m=A[k*nca+j]; 
            for(l=0; l<nca; l++){
                A[k*nca+l] = A[k*nca+l] - m*A[j*nca+l];
                C[k*nca+l] = C[k*nca+l] - m*C[j*nca+l];
            }
    }
    return 0;
}
/*============================================================================*/
int mta_2dswaprows(double *ar, const int r1, const int r2, const int nf, const int nc){
    if(nf<=0 || nc<=0) return -1;
    int i;
    double temp;
    for(i=0;i<nc;i++){
        temp = ar[r1*nc+i]; 
        ar[r1*nc+i] = ar[r2*nc+i]; 
        ar[r2*nc+i] = temp; 
    }
    return 0;
}
/*============================================================================*/
int mta_2dswapcols(double *ar, const int c1, const int c2, const int nf, const int nc){
    if(nf<=0 || nc<=0) return -1;
    int i;
    double temp;
    for(i=0;i<nf;i++){
        temp = ar[i*nc+c1]; 
        ar[i*nc+c1] = ar[i*nc+c2]; 
        ar[i*nc+c2] = temp; 
    }
    return 0;
}
/*============================================================================*/