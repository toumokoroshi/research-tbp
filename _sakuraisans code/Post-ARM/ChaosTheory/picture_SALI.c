#include <stdio.h>
#include <math.h>
#define epsilon 1.0e-10;
FILE* outputfile;
FILE* myfile;
double mu = 3.003e-6;
double C =  3.000100;
double y_min, y_max;
double k = -1;  ///////k=1:prograde motion, k=-1:retrograde motion
int i, j = 0;



int main() {
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
    fprintf(myfile, "set yrange[-0.01:0.01]\n");
    fprintf(myfile, "set cbrange[0:1.6]\n");
    fprintf(myfile, "set xlabel 'x-axis'\n");
    fprintf(myfile, "set ylabel 'y-axis'\n");
    fprintf(myfile, "set cblabel 'SALI'\n");
    fprintf(myfile, "set xlabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set ylabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set cblabel font 'Times New Roman, 20'\n");
    fprintf(myfile, "set pm3d map\n");
    fprintf(myfile, "set terminal png\n");
    fprintf(myfile, "set output 'k= -1.000000 ,C= 3.000100.png'\n");
    fprintf(myfile, "set palette defined (0.0 \"blue\", 0.1 \"green\", 0.2 \"yellow\",0.3 \"red\")\n");
    fprintf(myfile, "splot 'output_k= -1.000000 ,C= 3.000100.d' with pm3d\n");
    fprintf(myfile, "pause -1\n");
    return 0;
}
