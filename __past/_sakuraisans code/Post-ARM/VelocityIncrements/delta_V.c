#include<math.h>
#include<stdio.h>
#include<stdlib.h>
 
FILE* outputfile;
FILE* myfile;
///INPUT DATA
double a = 1.05;
double e = 0.05;
///Variable DATA
double cof, sif;
double x, y, V1, V2, C_pre, vr, vt;
double C;
double mu = 3.0e-6;
double delta_V, q1, q2;

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

int main() {
    outputfile = fopen("outputfile.d", "w");
    C_pre = 1 / a + 2 * sqrt(a * (1 - e * e));
    C = C_pre + 0.0005;
    if (outputfile == NULL) {
        printf("  can not open write file.");
    }
    printf("C_pre = %f, C = %f", C_pre, C);

    for (x = 0.97; x <= 1.03; x+=0.0001) {
        for (y = -0.03; y <= 0.03; y+=0.0001) {

            if ((x-1) * (x-1) + y * y > 0.03 * 0.03) {
                continue;
            }
            q1 = distance1(x, y);
            q2 = distance2(x, y);
            if (q2 < 0.0004) {
                fprintf(outputfile, "%f  %f \n", x, y);
                continue;
            }
            V1 = (x * x + y * y) + 2 * (1 - mu) / q1 + 2 * mu / q2 + mu * (1 - mu) - C_pre;
            V2 = (x * x + y * y) + 2 * (1 - mu) / q1 + 2 * mu / q2 + mu * (1 - mu) - C;
            if (V2 < 0) {
                fprintf(outputfile, "%f  %f  \n", x, y);
                continue;
            }
            V1 = sqrt(V1);
            V2 = sqrt(V2);
            delta_V = V1 - V2;
            fprintf(outputfile, "%f  %f  %f\n", x, y, delta_V);
        }
        fprintf(outputfile, "\n");
    }

    myfile = popen("gnuplot -persist", "w");
    fprintf(myfile, "unset key\n");
    fprintf(myfile, "set title 'a = %f,e = %f'\n", a, e);
    fprintf(myfile, "set size ratio 1 1\n");
    fprintf(myfile, "set xrange[0.97:1.03]\n");
    fprintf(myfile, "set yrange[-0.03:0.03]\n");
    //fprintf(myfile, "set cbrange[0:0.005]\n");
    fprintf(myfile, "set xlabel 'x axis'\n");
    fprintf(myfile, "set ylabel 'y axis'\n");
    fprintf(myfile, "set label 1 at graph 0.4,0.5 'L1'\n");
    fprintf(myfile, "set label 2 at graph 0.6,0.5 'L2'\n");
    fprintf(myfile, "set cblabel '{/Symbol D}V'\n");
    fprintf(myfile, "set pm3d map\n");
    fprintf(myfile, "set palette defined (0.0 \"blue\", 0.1 \"green\", 0.2 \"yellow\",0.3 \"red\")\n");
    fprintf(myfile, "splot 'outputfile.d' with pm3d\n");
    
    printf("\n Finish simulation!!!!!!\n");
    fclose(outputfile);
    return 0;
}