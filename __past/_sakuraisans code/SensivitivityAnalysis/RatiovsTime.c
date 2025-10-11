#include<math.h>
#include<stdio.h>
#include<stdlib.h>

FILE* outputfile;
FILE* myfile;
double Cj=2.999913;
double C_now;
double a=1.023;
double e=0.027;
double alfa;
double x=0.990200;
double y=0.000195;
double V, V1, delta_Vx, delta_Vy;
double vx, vy, vx1, vy1, C;
double q1,q2;
double X[4],Y[4],Y1[4],Y2[4],Y3[4],Z[4],K[4][6];
double mu = 3.003e-6;
double time;
double dt=0.0001;
double t_max = 300;
double k = -1;///k=1:posigrade motion, k=-1:retrograde motion
int i, j;
int count;
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


////equation of motion
// ////equation of motion
void equation_sy(double Y[4], const double c, const double d, double Z[4]){
       Z[0]=Y[0]+(+Y[1]+Y[2])*dt*c;
       Z[1]=Y[1]+(-Y[0]+Y[3])*dt*c;
       Z[2]=Y[2]+(+Y[3]-(1-mu)*(Z[0]+mu)/pow(((Z[0]+mu)*(Z[0]+mu)+Z[1]*Z[1]),3./2.)-mu*(Z[0]-1+mu)/pow(((Z[0]-1+mu)*(Z[0]-1+mu)+Z[1]*Z[1]),3./2.))*dt*d;
       Z[3]=Y[3]+(-Y[2]-(1-mu)*(Z[1])/pow(((Z[0]+mu)*(Z[0]+mu)+Z[1]*Z[1]),3./2.)-mu*(Z[1])/pow(((Z[0]-1+mu)*(Z[0]-1+mu)+Z[1]*Z[1]),3./2.))*dt*d;
}


int main() {
    C = 1/a+2.0*sqrt(a*(1-e*e));
    outputfile = fopen("outputfile_alfatime.d", "w");
    if (outputfile == NULL) {
        printf("  can not open write file.");
    }
    for(alfa=0.7 ; alfa <= 1.30 ; alfa += 0.01){

        // if (fabs(x - 1.0 + mu) < 0.00007) {///distance from the Earth
        //     fprintf(outputfile, "%f  %f  \n", x, C);
        //     continue;
        // }
        // else if ((x * x + 2*(1 - mu) / fabs(x + mu) + 2*mu / fabs(x - 1 + mu) + mu * (1 - mu)  - C) < 0) {///Out of ZVC
        //     fprintf(outputfile, "%f  %f  \n", x, C);
        //     continue;
        // }

        q1 = distance1(x, y);
        q2 = distance2(x, y);

        V = sqrt( (x * x + y * y) + 2 * (1 - mu) / q1 + 2 * mu / q2 + mu * (1 - mu) - Cj );
        V1 = sqrt( (x * x + y * y) + 2 * (1 - mu) / q1 + 2 * mu / q2 + mu * (1 - mu) - C );

        vx1 = -k * y /q2 * V1 ;
        vy1 = k * (x - 1 + mu) /q2 * V1 ;
        vx = -k * y /q2 * V ;
        vy = k * (x - 1 + mu) /q2 * V ;
        delta_Vx = vx - vx1;
        delta_Vy = vy - vy1;

        vx = vx1 + delta_Vx * alfa;
        vy = vy1 + delta_Vy * alfa;
        C_now = (x * x + y * y) + 2.0 * (1 - mu) / q1 + 2.0 * mu / q2 + mu * (1 - mu) - (vx * vx + vy * vy);
        time = 0.0;

        // ///initial condition for Symplectic method///
        X[0] = x;
        X[1] = y;
        X[2] = vx - y;
        X[3] = vy + x;
        /////End of initial setting for Symplectic method///

        while (time<t_max) {
            if (((X[0] - 1 + mu) * (X[0] - 1 + mu) + X[1] * X[1]) < 0.00007 * 0.00007) {///distance from the Earth50
                fprintf(outputfile, "%f  %f  %f\n", alfa, time, C_now);
                break;
            }else if (((X[0] - 1 + mu) * (X[0] - 1 + mu) + X[1] * X[1]) > 0.03 * 0.03) {///distance from the Earth0
                fprintf(outputfile, "%f  %f  %f\n", alfa, time, C_now);
                break;
            }
            ///////////////////////Runge-Kutta method////////////////////////////
            equation_sy(X,c1,d1,Y1);
            equation_sy(Y1,c2,d2,Y2);
            equation_sy(Y2,c3,d3,Y3);
            equation_sy(Y3,c4,d4,X);
            time += dt;

        }
        if (time >= t_max) {
            fprintf(outputfile, "%f  %f  %f\n", alfa, time, C_now);///100
        }
    }
    fclose(outputfile);
    printf("\nFinish simulation!!!!!\n");

    myfile = popen("gnuplot -persist", "w");
    fprintf(myfile, "unset key\n");
    if (k == 1) {
        fprintf(myfile, "set title 'Prograde motion'\n");
    }
    else if (k == -1) {
        fprintf(myfile, "set title 'Retrograde motion'\n");
    }
    fprintf(myfile, "set size ratio 1 1\n");
    fprintf(myfile, "set xrange[0.7:1.3]\n");
    //fprintf(myfile, "set yrange[2.995:3.015]\n");
    fprintf(myfile, "set xlabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set ylabel font 'Times New Roman, 20'\n");
    //fprintf(myfile, "set cblabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set xlabel 'Ratio, -'\n");
    fprintf(myfile, "set ylabel 'Time, -'\n");
    //fprintf(myfile, "set cblabel 'Time, -'\n");
    //fprintf(myfile, "set pm3d map\n");
    fprintf(myfile, "set terminal svg\n");
    fprintf(myfile, "set output 'Ratio-time.svg'\n",k);
    //fprintf(myfile, "set palette defined (0.0 \"blue\", 0.1 \"green\", 0.2 \"yellow\",0.3 \"red\")\n");
    fprintf(myfile, "plot 'outputfile_alfatime.d' w l\n");
    return 0;
}