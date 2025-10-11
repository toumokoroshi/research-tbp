#include<math.h>
#include<stdio.h>
#include<stdlib.h>

FILE* outputfile;
FILE* myfile;

double x, C;
double y, vy;
double X[4],Y[4],Y1[4],Y2[4],Y3[4],Z[4],K[4][6];
double mu = 3.003e-6;
double time;
double dt=0.0001;
double t_max = 10;
double x_min = 0.99;
double x_max = 1.01;
double x_step = 0.0001;
double k = 1;///k=1:posigrade motion, k=-1:retrograde motion
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


////equation of motion
// ////equation of motion
void equation_sy(double Y[4], const double c, const double d, double Z[4]){
       Z[0]=Y[0]+(+Y[1]+Y[2])*dt*c;
       Z[1]=Y[1]+(-Y[0]+Y[3])*dt*c;
       Z[2]=Y[2]+(+Y[3]-(1-mu)*(Z[0]+mu)/pow(((Z[0]+mu)*(Z[0]+mu)+Z[1]*Z[1]),3./2.)-mu*(Z[0]-1+mu)/pow(((Z[0]-1+mu)*(Z[0]-1+mu)+Z[1]*Z[1]),3./2.))*dt*d;
       Z[3]=Y[3]+(-Y[2]-(1-mu)*(Z[1])/pow(((Z[0]+mu)*(Z[0]+mu)+Z[1]*Z[1]),3./2.)-mu*(Z[1])/pow(((Z[0]-1+mu)*(Z[0]-1+mu)+Z[1]*Z[1]),3./2.))*dt*d;
}


int main() {

    outputfile = fopen("outputfile_xC_pro.d", "w");
    if (outputfile == NULL) {
        printf("  can not open write file.");
    }
    for(x=x_min ; x<=x_max ; x+=x_step){
        count++;
        printf("\b\b\b\b\b\b\b\b\b\b\b");
        printf("%d / %d", count, (int)((x_max - x_min) / x_step + 1));
        for (C = 2.995; C <= 3.015 ; C+=0.0001){
            if (fabs(x - 1.0 + mu) < 0.00007) {///distance from the Earth
                fprintf(outputfile, "%f  %f  \n", x, C);
                continue;
            }
            else if ((x * x + 2*(1 - mu) / fabs(x + mu) + 2*mu / fabs(x - 1 + mu) + mu * (1 - mu)  - C) < 0) {///Out of ZVC
                fprintf(outputfile, "%f  %f  \n", x, C);
                continue;
            }
            else if (x > (1 - mu)) {
                vy = k * sqrt(x * x + 2 * (1 - mu) / fabs(x + mu) + 2 * mu / fabs(x - 1 + mu) + mu * (1 - mu) - C);
            }
            else if (x < (1 - mu)) {
                vy = -k * sqrt(x * x + 2 * (1 - mu) / fabs(x + mu) + 2 * mu / fabs(x - 1 + mu) + mu * (1 - mu) - C);
            }
            time = 0.0;

            // ///initial condition for Symplectic method///
            X[0] = x;
            X[1] = 0.0;
            X[2] = 0.0;
            X[3] = vy+x;
            /////End of initial setting for Symplectic method///

            while (time<t_max) {
                if (((X[0] - 1 + mu) * (X[0] - 1 + mu) + X[1] * X[1]) < 0.00007 * 0.00007) {///distance from the Earth50
                    fprintf(outputfile, "%f  %f  %f\n", x, C, time);
                    break;
                }else if (((X[0] - 1 + mu) * (X[0] - 1 + mu) + X[1] * X[1]) > 0.03 * 0.03) {///distance from the Earth0
                    fprintf(outputfile, "%f  %f  %f\n", x, C, time);
                    break;
                }
                ///////////////////////Runge-Kutta method////////////////////////////
                equation_sy(X,c1,d1,Y1);
                equation_sy(Y1,c2,d2,Y2);
                equation_sy(Y2,c3,d3,Y3);
                equation_sy(Y3,c4,d4,X);
                time+=dt;

            }
            if (time >= t_max) {
                fprintf(outputfile, "%f  %f  %f\n", x, C, time);///100
            }
        }
        fprintf(outputfile, "\n");
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
    fprintf(myfile, "set xrange[0.99:1.01]\n");
    fprintf(myfile, "set yrange[2.995:3.015]\n");
    fprintf(myfile, "set cbrange[0:10]\n");
    fprintf(myfile, "set xlabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set ylabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set cblabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set xlabel 'x-axis'\n");
    fprintf(myfile, "set ylabel 'Jacobi integral C_{j}, -'\n");
    fprintf(myfile, "set cblabel 'Time, -'\n");
    fprintf(myfile, "set pm3d map\n");
    fprintf(myfile, "set terminal png\n");
    fprintf(myfile, "set output 'xC_k=%f.png'\n",k);
    fprintf(myfile, "set palette defined (0.0 \"blue\", 0.1 \"green\", 0.2 \"yellow\",0.3 \"red\")\n");
    fprintf(myfile, "splot 'outputfile_xC_pro.d' with pm3d\n");
    return 0;
}