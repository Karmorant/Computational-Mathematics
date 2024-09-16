#define M_PI acos(-1.0)

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>


union f64_union
{
        double d;
        uint64_t i;
} a, b ,c, d, e;

double FX_1(double x)
{
        double res;
        return res = sin(pow(x, 2)); /* cos(sin(x)); */ /* exp(sin(cos(x))); */  /* log(x + 3); */  /* pow(x + 3, 0.5) ; */
}


int main()
{

        FILE *fp_A = fopen("data_A.txt", "w");
        FILE *fp_B = fopen("data_B.txt", "w");
        FILE *fp_C = fopen("data_C.txt", "w");
        FILE *fp_D = fopen("data_D.txt", "w");
        FILE *fp_E = fopen("data_E.txt", "w");

        double x_0 = M_PI/6;
        double fx_0 =  2 * cos(pow(x_0, 2)) * x_0;  /* -sin(sin(x_0)) * cos(x_0); */ /*  -exp(sin(cos(x_0))) * cos(cos(x_0)) * sin(x_0); */ /* pow(x_0 + 3, -1); */ /*  0.5 * pow(x_0 + 3, -0.5); */

        a.i = 0x3FF5555555555555; // 4/3
        b.i = 0x3FD5555555555555; // 1/3
        c.i = 0x3FF8000000000000; // 3/2
        d.i = 0x3FE3333333333333; // 3/5
        e.i = 0x3FB999999999999A; // 1/10

        printf("%f\n", FX_1(x_0 + 1));
        printf("%f\n", FX_1(1));
        printf("%f\n", (FX_1(1 + 1) - FX_1(1)) / 1);
        printf("%f\n", fx_0);
        printf("%f\n", ((a.d) * ((FX_1(x_0 + 1) - FX_1(x_0 - 1)) / (2 * 1))));

        for (int x = 1; x <= 21; ++x)
        {
                double h = pow(2, 1 - x);
                double fa = (FX_1(x_0 + h) - FX_1(x_0)) / h;
                double fb = (FX_1(x_0) - FX_1(x_0 - h)) / h;
                double fc = (FX_1(x_0 + h) - FX_1(x_0 - h)) / (2 * h);
                double fd = ((a.d) * ((FX_1(x_0 + h) - FX_1(x_0 - h)) / (2 * h))) - ((b.d) * ((FX_1(x_0 + 2 * h) - FX_1(x_0 - 2 * h)) / (4 * h)));
                double fe = ((c.d) * ((FX_1(x_0 + h) - FX_1(x_0 - h)) / (2 * h))) - ((d.d) * ((FX_1(x_0 + 2 * h) - FX_1(x_0 - 2 * h)) / (4 * h))) + ((e.d) * ((FX_1(x_0 + 3 * h) - FX_1(x_0 - 3 * h)) / (6 * h)));

                fprintf(fp_A, "%f %.23f\n", h, fabs(fx_0 - fa));
                fprintf(fp_B, "%f %.23f\n", h, fabs(fx_0 - fb));
                fprintf(fp_C, "%f %.23f\n", h, fabs(fx_0 - fc));
                fprintf(fp_D, "%f %.23f\n", h, fabs(fx_0 - fd));
                fprintf(fp_E, "%f %.23f\n", h, fabs(fx_0 - fe));
        }

        fclose(fp_A);
        fclose(fp_B);
        fclose(fp_C);
        fclose(fp_D);
        fclose(fp_E);

        FILE *gnuplotPipe = popen("./plot.pgi", "w");
        fclose(gnuplotPipe);

        return 0;
}