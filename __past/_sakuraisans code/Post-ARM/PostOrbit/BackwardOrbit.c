#include<math.h>
#include<stdio.h>
#include<stdlib.h>
 
FILE* outputfile;
FILE* myfile;
FILE* myfile2;
///INPUT DATA
double a = 1.062;
double e = 0.052;
///Variable DATA
double Y[4], Z[4], Z1[4], K[4][6], Y1[4],Y2[4],Y3[4];
double X[4];
double mu = 3.003e-6;
double t_end = 50;
double dt = 0.0001;
double t = 0.0;
double ecc, C_pre;
double cof, sif;
double x, y, V1;
double delta_V, q1, q2;
double q1, q2;
double x, y, vy, vx;
double V;
double k = -1;  ///////k=1:prograde motion, k=-1:retrograde motion
int count;
int i, j = 0;
double c1=1/(2*(2-pow(2,1/3)));
double c2=(1-pow(2,1/3))/(2*(2-pow(2,1/3)));
double c3=(1-pow(2,1/3))/(2*(2-pow(2,1/3)));
double c4=1/(2*(2-pow(2,1/3)));
double d1=1/(2-pow(2,1/3));
double d2=-pow(2,1/3)/(2-pow(2,1/3));
double d3=1/(2-pow(2,1/3));
double d4=0.0;



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

// ////equation of motion of backward 
void equation_sy(double Y[4], const double c, const double d, double Z[4]){
       Z[0]=Y[0]-(+Y[1]+Y[2])*dt*c;
       Z[1]=Y[1]-(-Y[0]+Y[3])*dt*c;
       Z[2]=Y[2]-(+Y[3]-(1-mu)*(Z[0]+mu)/pow(((Z[0]+mu)*(Z[0]+mu)+Z[1]*Z[1]),3./2.)-mu*(Z[0]-1+mu)/pow(((Z[0]-1+mu)*(Z[0]-1+mu)+Z[1]*Z[1]),3./2.))*dt*d;
       Z[3]=Y[3]-(-Y[2]-(1-mu)*(Z[1])/pow(((Z[0]+mu)*(Z[0]+mu)+Z[1]*Z[1]),3./2.)-mu*(Z[1])/pow(((Z[0]-1+mu)*(Z[0]-1+mu)+Z[1]*Z[1]),3./2.))*dt*d;
}


int main(){
    printf("\n  **********************************************************************\n");
    printf("  *   Calculator of Orbit before mission with Symplectic method        *\n");
    printf("  *                                                                    *\n");
    printf("  *   This source code was written by Yuto Sakurai,                    *\n");
    printf("  *         Department of Aerospace Engineering,  Nagoya University    *\n");
    printf("  *                                                                    *\n");
    printf("  *  Last Update : 2020.10.18                                          *\n");
    printf("  **********************************************************************\n\n");

    outputfile = fopen("output_backwardorbit.d", "w");
    if (outputfile == NULL) {
        printf("  Can not open write file");
    }
    printf(" Simulating  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");

    C_pre = 1 / a + 2 * sqrt(a * (1 - e * e));

    ///input capture position////
    x=0.99020;
    y=0.000195;

    q1 = distance1(x, y);
    q2 = distance2(x, y);

    V = sqrt(x * x + y * y + 2 * (1 - mu) / q1 + 2 * mu / q2 + mu * (1 - mu) - C_pre);
    vx = -k * V * y / q2;
    vy = k * V * (x - 1 + mu) / q2;


    // ///initial condition for Symplectic method///
    t = 0.0;

    X[0] = x;
    X[1] = y;
    X[2] = vx-y;
    X[3] = vy+x;

    /////End of initial setting for Symplectic method///

    while (t < t_end) {
        t += dt;
        q2 = distance2(X[0], X[1]);
        if (q2 < 0.00007) {
            break;
        }
        // }else if(q2 > 0.03){
        //     break;
        // }

        /////////////////////////////////////////////////Symplectic Integrator/////////////////
        equation_sy(X,c1,d1,Y1);
        equation_sy(Y1,c2,d2,Y2);
        equation_sy(Y2,c3,d3,Y3);
        equation_sy(Y3,c4,d4,X);

        q2 = distance2(X[0], X[1]);
        V1=(X[2]+X[1])*(X[2]+X[1])+(X[3]-X[0])*(X[3]-X[0]);
        V1 = sqrt(V1);
        V = V1 + k * q2;
        ecc = q2 * V * V / mu - 1;

        fprintf(outputfile, "%f  %f  %f  %f\n",t, ecc, X[0], X[1]);
        ///////////////////////////////////////////////End Symplectic Integrator///////////////////
    }
    printf("\n Finish simulation!!!!!!\n");
    fclose(outputfile);



    myfile = popen("gnuplot -persist", "w");
    //fprintf(myfile, "unset key\n");
    if (k == 1) {
        fprintf(myfile, "set title 'Prograde motion'\n");
    }
    else if (k == -1) {
        fprintf(myfile, "set title 'Retrograde motion'\n");
    }
    fprintf(myfile, "set size ratio 1 1\n");
    // fprintf(myfile, "set xrange[0.97:1.03]\n");
    // fprintf(myfile, "set yrange[-0.03:0.03]\n");
    fprintf(myfile, "set xlabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set ylabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set label 1 point pt 7 ps 2 lc rgb 'blue' at 1 ,0\n");
    fprintf(myfile, "set xlabel 'x axis'\n");
    fprintf(myfile, "set ylabel 'y axis'\n");
    fprintf(myfile, "set terminal png\n");
    fprintf(myfile, "set output 'Orbit.png'\n");
    fprintf(myfile, "plot 'output.d' using 3:4 w l title ' Captured Orbit', 'output_backwardorbit.d' using 3:4 w l title 'Natural Orbit',  'output.d' using 3:4 every ::0::0  pt 7 ps 2 title 'Capture point'\n");
    //fprintf(myfile, "pause -1\n");

    myfile2 = popen("gnuplot -persist", "w");
    fprintf(myfile2, "unset key\n");
    if (k == 1) {
        fprintf(myfile2, "set title 'Prograde motion'\n");
    }
    else if (k == -1) {
        fprintf(myfile2, "set title 'Retrograde motion'\n");
    }
    fprintf(myfile2, "set size ratio 1 1\n");
    fprintf(myfile2, "set xrange[0:20]\n");
    //fprintf(myfile2, "set yrange[0:1.2]\n");
    fprintf(myfile2, "set xlabel font 'Times New Roman, 20'\n");
    fprintf(myfile2, "set ylabel font 'Times New Roman, 20'\n");
    fprintf(myfile2, "set xlabel 'time, -'\n");
    fprintf(myfile2, "set ylabel 'eccentricity, -'\n");
    fprintf(myfile2, "set terminal png\n");
    fprintf(myfile2, "set output 'Pre_eccentricity.png'\n");
    fprintf(myfile2, "plot 'output_backwardorbit.d' using 1:2 w l\n");
    //fprintf(myfile2, "pause -1\n");

    return 0;
}