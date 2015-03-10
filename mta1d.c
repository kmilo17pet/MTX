#include "mta1d.h"

/*============================================================================*/
void mta_update_reg(double regarray[],const double newp,const int length){
        int i;
        for (i=length-1;i>=1;i=i-1)
            regarray[i]=regarray[i-1];
        regarray[0]=newp;
}
/*============================================================================*/
double mta_median_filt(double regarray[],const double newp,const int length){
        int i;
        double accum=0;
        for (i=length-1;i>=1;i=i-1)
        {
            regarray[i]=regarray[i-1];
            accum+=regarray[i];
        }
        regarray[0]=newp;
        return (accum+newp)/length;
}
/*============================================================================*/
void mta_1dinit(double ar[],const int arl,const double val){
    int k;
    for(k=0;k<arl;k++) ar[k]=val;
}
/*============================================================================*/
double mta_1dpolyval(double ar[],const int arl,const double val){
    int k;
    double outv=0;
    for (k=0;k<arl;k++) outv+=ar[arl-k-1]*pow(val,k);
    return outv;
}
/*============================================================================*/
double mta_1dsum(const double ar[],const int arl){
    int k;
    double asum=0;
    for (k=0;k<arl;k++) asum+=ar[k];
    return asum;
}
/*============================================================================*/
double mta_1dprod(const double ar[],const int arl){
    int k;
    double aprd=1;
    for (k=0;k<arl;k++) aprd*=ar[k];
    return aprd;
}
/*============================================================================*/
int mta_1doper(double out[],const double a[],const double b[],const int lout,const int la,const int lb,const int oper){
    if (((la!=lb)||(lout!=la))||(lout!=lb)) return -1;
    int numel=(la<lb)? la:lb;
    numel=(numel<lout)? numel:lout;
    int k;
    switch (oper){
        case 1:   for(k=0;k<numel;k++) out[k]=a[k]+b[k];
        break;
        case 2:   for(k=0;k<numel;k++) out[k]=a[k]-b[k];
        break;
        case 3:   for(k=0;k<numel;k++) out[k]=a[k]*b[k];
        break;
        case 4:   for(k=0;k<numel;k++) out[k]=a[k]/b[k];
        break;
        default:  for(k=0;k<numel;k++) out[k]=a[k]+b[k];
        break;
    }
}
/*============================================================================*/
void mta_1dminus_one(double ar[],const int l){
    int k;
    for (k=0;k<l;k++) ar[k]=-ar[k];
}
/*============================================================================*/
