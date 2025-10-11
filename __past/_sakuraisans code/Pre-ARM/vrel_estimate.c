#include<math.h>
#include<stdio.h>
#include<stdlib.h>
 
    FILE* outputfile;
    FILE* outputfile1;
    FILE* myfile;
    FILE* myfile1;
    FILE* myfile2;
    double V1, vr, vt, angle;
    double a, e;
    double cof, sif;
    double t20 = tan(20 * 3.141592 / 180);
    double t40 = tan(40 * 3.141592 / 180);
    double t60 = tan(60 * 3.141592 / 180);
    double t80 = tan(80 * 3.141592 / 180);

    int main() {
        
        outputfile = fopen("outputfile.d", "w");
        outputfile1 = fopen("outputfile1.d", "w");
        if (outputfile == NULL) {
            printf("  can not open write file.");
        }
        for (a = 0; a <= 3; a+=0.001) {
            for (e = 0; e <= 1; e+=0.001) {
                

                cof = (a * (1 - e * e) - 1) / e;
                if (cof * cof > 1) {
                    //fprintf(outputfile, "%f  %f  0.0  \n", a, e);
                    continue;
                }
                sif = sqrt(1 - cof * cof);

                vr = e * sif / sqrt(a * (1 - e * e));
                vt = (1 + e * cof) / sqrt(a * (1 - e * e));

                V1 = (vt - 1) * (vt - 1) + vr * vr;
                V1 = sqrt(V1);
                angle = atan(vt / vr);

                fprintf(outputfile, "%f  %f  %f\n", a, e, V1);
            }
            fprintf(outputfile, "\n");
        }
        fclose(outputfile);

        for (e = 0; e <= 1.0; e += 0.001) {
            fprintf(outputfile1, "%f  %f", e, ((2 * t20 * t20) + sqrt((2 * t20 * t20) * (2 * t20 * t20) - 4 * t20 * t20 * (1 - e * e) * (t20 * t20 * +1))) / (2 * (1 - e * e) * (t20 * t20 + 1)));
            fprintf(outputfile1, "  %f", ((2 * t40 * t40) + sqrt((2 * t40 * t40) * (2 * t40 * t40) - 4 * t40 * t40 * (1 - e * e) * (t40 * t40 * +1))) / (2 * (1 - e * e) * (t40 * t40 + 1)));
            fprintf(outputfile1, "  %f", ((2 * t60 * t60) + sqrt((2 * t60 * t60) * (2 * t60 * t60) - 4 * t60 * t60 * (1 - e * e) * (t60 * t60 * +1))) / (2 * (1 - e * e) * (t60 * t60 + 1)));
            fprintf(outputfile1, "  %f", ((2 * t80 * t80) + sqrt((2 * t80 * t80) * (2 * t80 * t80) - 4 * t80 * t80 * (1 - e * e) * (t80 * t80 * +1))) / (2 * (1 - e * e) * (t80 * t80 + 1)));
            fprintf(outputfile1, "\n");
        }

        myfile2 = popen("gnuplot -persist", "w");
        fprintf(myfile2, "set size ratio 1 1\n");
        fprintf(myfile2, "set xrange[0:3]\n");
        fprintf(myfile2, "set yrange[0:1]\n");
        fprintf(myfile2, "set xlabel 'Semi-major axis a_{ast}(au)'\n");
        fprintf(myfile2, "set ylabel 'Eccentricity e_{ast}'\n");
        fprintf(myfile2, "plot 'outputfile1.d' using 2:1 w l title '20 deg','outputfile1.d' using 3:1 w l title '40 deg','outputfile1.d' using 4:1 w l title '60 deg','outputfile1.d' using 5:1 w l title '80 deg'\n");


        myfile = popen("gnuplot -persist", "w");        
        fprintf(myfile, "set size ratio 1 1\n");
        fprintf(myfile, "set xrange[0:3]\n");
        fprintf(myfile, "set yrange[0:1]\n");
        fprintf(myfile, "set zrange[0:2]\n");
        fprintf(myfile, "set key at 2.85,0.42\n");
        fprintf(myfile, "set key box linetype 1 linewidth 1 lc rgb 'black'\n");
        fprintf(myfile, "set label 1 at graph 0.1,0.25 'No orbital intersection' rotate by 90\n");
        fprintf(myfile, "set label 2 at graph 0.48,0.10 'No orbital intersection' rotate by 0\n");
        //fprintf(myfile, "set label 3 at graph 0.35,0.9 '20 deg' rotate by 15\n");
        //fprintf(myfile, "set label 4 at graph 0.35,0.6 '40 deg' rotate by 43\n");
        //fprintf(myfile, "set label 5 at graph 0.35,0.3 '60 deg' rotate by 70\n");
        //fprintf(myfile, "set label 6 at graph 0.35,0.15 '80 deg' rotate by 80\n");
        fprintf(myfile, "set xlabel font 'Times New Roman, 20'\n");
        fprintf(myfile, "set ylabel font 'Times New Roman, 20'\n");
        fprintf(myfile, "set cblabel font 'Times New Roman, 20'\n");
        fprintf(myfile, "set xlabel 'Semi-major axis, au'\n");//{/Arial-Italic a}_{ast},
        fprintf(myfile, "set ylabel 'Eccentricity, -'\n");
        fprintf(myfile, "set cblabel 'Relative velocity, -'\n");
        fprintf(myfile, "set pm3d map\n");
        //fprintf(myfile, "set palette defined (0.0 \"blue\", 0.3 \"green\", 0.6 \"yellow\", 0.8 \"red\")\n");
        fprintf(myfile, "splot 'outputfile.d' with pm3d title '', 'outputfile1.d' using 2:1:($2-$2) with lines dt '.' lw 3 title '20 deg' , 'outputfile1.d' using 3:1:($2-$2) with lines dt '-' lw 3 title '40 deg' ,'outputfile1.d' using 4:1:($2-$2) with lines dt '_' lw 3 title '60 deg','outputfile1.d' using 5:1::($2-$2) with lines title '80 deg' lw 3\n");

//, 'outputfile1.d' using 2:1:($2-$2) with lines dt '.' lw 3 title '20 deg' , 'outputfile1.d' using 3:1:($2-$2) with lines dt '-' lw 3 title '40 deg' ,'outputfile1.d' using 4:1:($2-$2) with lines dt '_' lw 3 title '60 deg','outputfile1.d' using 5:1::($2-$2) with lines title '80 deg' lw 3

        myfile1 = popen("gnuplot -persist", "w");
        fprintf(myfile1, "unset surface\n");
        fprintf(myfile1, "set contour\n");
        fprintf(myfile1, "set view 0,0,1,1\n");
        fprintf(myfile1, "set size ratio 1 1\n");
        fprintf(myfile1, "set xrange[0:3]\n");
        fprintf(myfile1, "set yrange[0:1]\n");
        fprintf(myfile1, "set xlabel 'Semi-major axis (au)'\n");
        fprintf(myfile1, "set ylabel 'Eccentricity'\n");
        fprintf(myfile1, "set isosamples 100\n");
        fprintf(myfile1, "splot atan(x*(1-y*y)/sqrt(y*y-(x*(1-y*y)-1)))*180/3.14159\n");
       
        printf("\n Finish simulation!!!!!!\n");
        
        pclose(myfile);
        pclose(myfile1);
        pclose(myfile2);
        return 0;
    }