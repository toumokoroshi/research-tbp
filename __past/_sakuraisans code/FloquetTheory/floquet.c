#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define     EPS   1.0e-1     //  許容誤差
#define     N     4          //  行列の元数

FILE* outputfile;
FILE* myfile;
FILE* myfile2;
///Variable DATA
double X[4][4], next_X[4][4], Y[4], next_Y[4];
double mu = 3.003e-6;
double t_end = 3.518000;  ////1period time for monodromy matrix
double dt = 0.001;
double t = 0.0;
double C = 2.999913;////Input Jacobi integral//
double q1, q2;
double x, y, vy, vx;
int count;
double A[4][4];
void disp_eigenvalue(double a[4][4]);
void S_Householder(double *,int);
void S_QRmethod(double *,int);
void S_OutMat( double *, int, int );


double distance1(double x, double y) {
    q1 = (x + mu) * (x + mu) + y * y;
    q1 = sqrt(q1);
    return q1;
}

double distance2(double x, double y) {
    q2 = (x - 1 + mu) * (x - 1 + mu) + y * y;
    q2 = sqrt(q2);
    return q2;
}


void Equation(const double X[4][4], const double A[4][4], double K[4][4])
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++){K[i][j] = 0.0;}
    }
    
    {
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                for(int k = 0; k < 4; k++) { K[i][j] += A[i][k] * X[k][j]; }
            }
        }
    }
}

void Dynamics(const double Y[4], double K[4])
{
    K[0] = Y[2];
    K[1] = Y[3];
    K[2] = 2 * Y[3] + Y[0] - (1 - mu) * (Y[0] + mu) / pow(((Y[0] + mu) * (Y[0] + mu) + Y[1] * Y[1]), 3. / 2.) - mu * (Y[0] - 1 + mu) / pow(((Y[0] - 1 + mu) * (Y[0] - 1 + mu) + Y[1] * Y[1]), 3. / 2.);
    K[3] = -2 * Y[2] + Y[1] - (1 - mu) * (Y[1]) / pow(((Y[0] + mu) * (Y[0] + mu) + Y[1] * Y[1]), 3. / 2.) - mu * (Y[1]) / pow(((Y[0] - 1 + mu) * (Y[0] - 1 + mu) + Y[1] * Y[1]), 3. / 2.);
}


void RungeOneStep(const double x[4][4], const double A[4][4], double next_x[4][4])
{
    double k1[4][4], k2[4][4], k3[4][4], k4[4][4];
    double xk2[4][4], xk3[4][4], xk4[4][4];

    Equation(x, A, k1);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++){xk2[i][j] = x[i][j] + (dt / 2.0) * k1[i][j];}
    }
    Equation(xk2, A, k2);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++){xk3[i][j] = x[i][j] + (dt / 2.0) * k2[i][j];}
    }
    Equation(xk3, A, k3);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++){xk4[i][j] = x[i][j] + dt * k3[i][j];}
    }
    Equation(xk4, A, k4);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++){next_x[i][j] = x[i][j] + (dt / 6.0)*(k1[i][j] + 2.0 * k2[i][j] + 2.0 * k3[i][j] + k4[i][j]);}
    }
}


void RungeOneStep4Motion(double y[4], double next_y[4])
{
    double k1[4], k2[4], k3[4], k4[4];
    double yk2[4], yk3[4], yk4[4];

    Dynamics(y, k1);
    for(int i = 0; i < 4; i++){yk2[i] = y[i] + (dt / 2.0) * k1[i];}
    Dynamics(yk2, k2);
    for(int i = 0; i < 4; i++){yk3[i] = y[i] + (dt / 2.0) * k2[i];}
    Dynamics(yk3, k3);
    for(int i = 0; i < 4; i++){yk4[i] = y[i] + dt * k3[i];}
    Dynamics(yk4, k4);
    for(int i = 0; i < 4; i++){next_y[i] = y[i] + (dt / 6.0)*(k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);}
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
    int count=0;
    //   収束判定
    while( m != 1 ){
        
        count++;
        printf("\b\b\b\b\b\b\b\b");
        printf("%d",count);
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


int main() {
    printf("\n  **********************************************************************\n");
    printf("  *   Calculator of monodromy matrix with RK4                          *\n");
    printf("  *                                                                    *\n");
    printf("  *   This source code was written by Yuto Sakurai,                    *\n");
    printf("  *         Department of Aerospace Engineering,  Nagoya University    *\n");
    printf("  *                                                                    *\n");
    printf("  *   Last Update : 2020.12.16                                         *\n");
    printf("  **********************************************************************\n\n");

    ///input initial condition////
    Y[0] = 0.990200;
    Y[1] = 0.000198;
    Y[2] = 0.000630;
    Y[3] = 0.031311;

    outputfile = fopen("output.d", "w");
    if (outputfile == NULL) {
        printf("  Can not open write file");
    }
    t = 0.0;

    ///Xの初期値設定
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4 ; j++){ X[i][j] = (i == j)?1.0:0.0; }
    }

    while (t < t_end) {

        x = Y[0];
        y = Y[1];
        vx = Y[2];
        vy = Y[3];
        q1 = distance1(x, y);
        q2 = distance2(x, y);

        ///Calculate A matrix///
        A[0][0] = 0.0;     
        A[0][1] = 0.0;     
        A[0][2] = 1.0;      
        A[0][3] = 0.0;

        A[1][0] = 0.0;     
        A[1][1] = 0.0;     
        A[1][2] = 0.0;      
        A[1][3] = 1.0;

        A[2][0] = 1.0 - (1-mu) * (pow(q1, 2) - 3.0 * (x+mu) * (x+mu))/pow(q1, 5) - mu * (pow(q2, 2) - 3.0 * (x+mu-1) * (x+mu-1))/pow(q2, 5);//     
        A[2][1] = - (1-mu) * (- 3.0 * (x+mu) * y)/pow(q1, 5) - mu * (- 3.0 * (x+mu-1) * y)/pow(q2, 5);//      
        A[2][2] = 0.0;      
        A[2][3] = 2.0;

        A[3][0] = - (1-mu) * (- 3.0 * (x+mu) * y)/pow(q1, 5) - mu * (- 3.0 * (x+mu-1) * y)/pow(q2, 5);//     
        A[3][1] = 1.0 - (1-mu) * (pow(q1, 2) - 3.0 * y * y)/pow(q1, 5) - mu * (pow(q2, 2) - 3.0 * y * y)/pow(q2, 5);//     
        A[3][2] = -2.0;      
        A[3][3] = 0.0;       

        RungeOneStep(X, A, next_X);
        t += dt;
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++){ X[i][j] = next_X[i][j];}
        }

        RungeOneStep4Motion(Y, next_Y);
        for(int i = 0; i < 4; i++){ Y[i] = next_Y[i];}

        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                printf("%f      ", X[i][j]);
            }
            printf("\n");
        }
        printf("\n\n");

    }
    // for(int i = 0; i < 4; i++)
    // {
    //     for(int j = 0; j < 4; j++)
    //     {
    //         printf("%f      ", X[i][j]);
    //     }
    //     printf("\n");
    // }


    double a[N*N]={X[0][0],X[0][1],X[0][2],X[0][3], X[1][0],X[1][1],X[1][2],X[1][3], X[2][0],X[2][1],X[2][2],X[2][3], X[3][0],X[3][1],X[3][2],X[3][3]};
    //double a[N*N]={0.2,.32,.12,.3, .1,.15,.24,.32, .2,.24,.46,.36, .6,.4,.32,.2};
    //double a[N*N]={1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    printf(" Monodromy Matrix \n");
    S_OutMat( a, N, N );
	
    S_Householder( a, N );
    printf(" Matrix after Householder Transformation \n");
    S_OutMat( a, N, N );

    printf(" Matrix after QR Transformation \n");
    S_QRmethod( a, N );
    S_OutMat( a, N, N);


    printf("\n Finish simulation!!!!!!\n\n");
    fclose(outputfile);

    return 0;
}
