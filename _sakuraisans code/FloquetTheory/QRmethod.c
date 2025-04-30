
/**********************************************************************/
/*    QR法による固有値解析                                            */
/*            ファイル名：QRmethod.c                                  */
/******************************* programed by K.Minemura **************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define     EPS   1.0e-8     //  許容誤差
#define     N     4          //  行列の元数

void S_Householder(double *,int);
void S_QRmethod(double*,int);
void S_OutMat( double *, int, int );

int main( void ){
    //  ===== 行列(入力データ) ================================================
    //   double a[N*N]={4.0,-6.0,5.0, -6.0,3.0,4.0, 5.0,4.0,-3.0}; 
    //   double a[N*N]={-1.,2.,3.,3., 2.,-3.,4.,1., 3.,4.,-1.,2., 3.,1.,2.,-3.};
    double a[N*N]={0.2,.32,.12,.3, .1,.15,.24,.32, .2,.24,.46,.36, .6,.4,.32,.2};
    //   =======================================================================

    printf(" 行列 \n");
    S_OutMat( a, N, N );
	
    S_Householder( a, N );
    printf(" ハウスホルダー変換後の行列 \n");
    S_OutMat( a, N, N );

    printf(" QR法変換後の行列（対角項が固有値） \n");
    S_QRmethod( a, N );
    S_OutMat( a, N, N);

    return 0;
}
/******** Householder変換による上Hessenberg行列への返還副関数 **************/
/*    入出力:  a[n×n]=係数行列，n=正方行列の元数                          */
/***************************************************************************/
void S_Householder(double a[],int n){
    int i,j,k;
    double sum,sigma,v_norm,ud,uds,*u,*d,*ds;
    u = (double *)calloc(n, sizeof(double));
    d = (double *)calloc(n, sizeof(double));
    ds =(double *)calloc(n, sizeof(double));

    for(k=0; k<=n-3; k++){
       for(i=0; i<=k; i++)     u[i]=0.0;
       for(i=k+1; i<n; i++)    u[i]=a[n*i+k];

       //  変換行列 H の構築
       sum = 0.0;
       for(i=k+1; i<n; i++) sum=sum+u[i]*u[i];
       if( fabs(u[k+1]) < EPS ) continue;
       sigma = sqrt(sum)*u[k+1]/fabs(u[k+1]);
       u[k+1] += sigma;
       v_norm = sqrt(2.0*sigma*u[k+1]);
       for(i=k+1; i<n; i++) u[i] /= v_norm;

       //  相似変換
       for(i=0; i<n; i++){
          d[i] = 0.0; ds[i] = 0.0;
          for(j=k+1; j<=n-1; j++){
             d[i]  += a[n*i+j]*u[j];
             ds[i] += a[n*j+i]*u[j];
          }
       }
       ud = 0.0;  uds = 0.0;
       for(i=k+1; i<n; i++){
          ud  += u[i]*d[i];
          uds += u[i]*ds[i];
       }
       for(i=0; i<n; i++){
          d[i] = 2.0*(d[i]- ud*u[i]);
          ds[i]= 2.0*(ds[i]-uds*u[i]);
       }
       for(i=0; i<n; i++){
          for(j=0; j<n;j++)  a[n*i+j] -= u[i]*ds[j]+d[i]*u[j];
       }
    }
}
/********** QR分解による固有値解析する副関数 **************************/
/*    入出力: a[n*n]=正方行列     n=正方行列の元数                    */
/**********************************************************************/
void S_QRmethod(double a[],int n){
    int i,j,k,m;
    double a00,a01,a10,a11,lam1,lam2,sum1,sum2,wa,mu,sinx,cosx,*q,*w;
    q=(double *)calloc(n*n, sizeof(double));
    w=(double *)calloc(n,   sizeof(double));

    m = n;
    //   収束判定
    while( m != 1 ){
       if(fabs(a[n*(m-1)+m-2]) < EPS){ m=m-1; continue;}

       //   原点移動  mu
       a00 = a[n*(m-2)+m-2];   a01 = a[n*(m-2)+m-1];
       a10 = a[n*(m-1)+m-2];   a11 = a[n*(m-1)+m-1];
       sum1 = a00+a11;         sum2 = a00*a11-a01*a10;
       wa  = sum1*sum1-4.0*sum2;
       if(wa < 0.0) wa=0.0; else wa=sqrt(wa);
       lam1 = 0.5*(sum1+wa);   lam2=sum2/lam1;
       if(fabs(a11-lam1) < fabs(a11-lam2))  mu = a11-lam1;
       else  mu = a11-lam2;
       for(i=0; i<m; i++)  a[n*i+i] -= mu;

       //   QR分解
       for(i=0; i<m*m; i++)  q[i] = 0.0;
       for(i=0; i<m; i++)    q[m*i+i] = 1.0;
       for(i=0; i<m-1; i++){
          sum1 = sqrt(a[n*i+i]*a[n*i+i]+a[n*i+n+i]*a[n*i+n+i]);
          if(fabs(sum1) < EPS){
             sinx = 0.0;  cosx = 0.0;
          }else{
             sinx = a[n*i+n+i]/sum1;  cosx = a[n*i+i]/sum1;
          }
          for(j=i+1; j<m; j++){
             sum2 = a[n*i+j]*cosx+a[n*i+n+j]*sinx;
             a[n*i+n+j] = -a[n*i+j]*sinx+a[n*i+n+j]*cosx;
             a[n*i+j] = sum2;
          }
          a[n*i+n+i] = 0.0;
          a[n*i+i] = sum1;
          for(j=0; j<m; j++){
            sum2 = q[m*j+i]*cosx+q[m*j+i+1]*sinx;
             q[m*j+i+1] = -q[m*j+i]*sinx+q[m*j+i+1]*cosx;
             q[m*j+i] = sum2;
          }
       }
       for(i=0; i<m; i++){
          for(j=i; j<m; j++) w[j]=a[n*i+j];
          for(j=0; j<m; j++){
             sum1 = 0.0;
             for(k=i; k<m; k++) sum1 += w[k]*q[m*k+j];
             a[n*i+j] = sum1;
          }
       }
       for(i=0; i<m; i++) a[n*i+i] += mu;
    }
}
/*********  行列を表示する副関数  *******************************/
/*       入力；   a= 行列, n= 表示する列数, m= 表示する行数     */
/****************************************************************/
void S_OutMat( double a[], int n, int m )
{
    int i,j;
    for(i=0;i<n;i++){
        for(j=0; j<m; j++)  printf("  %10.6f",a[i*m+j]);
        printf("\n");
    }
}