#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BETA 10.0
#define ALPHA 1000.0


void f(double t, double *y, double *fval) {
    (void)t;
    fval[0] = 1.0 - y[0] * y[1];
    fval[1] = ALPHA * y[1] * (y[0] - (1.0 + BETA)/(y[1] + BETA));
}

void rosenbrock4(double t, double *y, double h, int n, double *ynew) {
    double k1[n], k2[n], k3[n], k4[n];
    double ytemp[n], fval[n];
    
    const double gamma = 0.5728160624821348;
    const double a21 = 1.0;
    const double a31 = 0.5;
    const double a32 = 0.5;
    const double a41 = 0.5;
    const double a42 = 0.5;
    const double a43 = 1.0;
    
    const double c21 = -1.0;
    const double c31 = -0.5;
    const double c32 = -0.5;
    const double c41 = -0.5;
    const double c42 = -0.5;
    const double c43 = -1.0;
    
    const double b1 = 1.0/6.0;
    const double b2 = 1.0/3.0;
    const double b3 = 1.0/3.0;
    const double b4 = 1.0/6.0;
    

    f(t, y, fval);
    for (int i = 0; i < n; i++) {
        k1[i] = fval[i] / (1.0 - gamma*h*c21);
        ytemp[i] = y[i] + a21*k1[i]*h;
    }
    

    f(t + 0.5*h, ytemp, fval);
    for (int i = 0; i < n; i++) {
        k2[i] = (fval[i] + h*c31*k1[i]) / (1.0 - gamma*h*c32);
        ytemp[i] = y[i] + (a31*k1[i] + a32*k2[i])*h;
    }
    

    f(t + 0.5*h, ytemp, fval);
    for (int i = 0; i < n; i++) {
        k3[i] = (fval[i] + h*(c41*k1[i] + c42*k2[i])) / (1.0 - gamma*h*c43);
        ytemp[i] = y[i] + (a41*k1[i] + a42*k2[i] + a43*k3[i])*h;
    }
    

    f(t + h, ytemp, fval);
    for (int i = 0; i < n; i++) {
        k4[i] = (fval[i] + h*(c41*k1[i] + c42*k2[i] + c43*k3[i])) / (1.0 - gamma*h*c43);
        ynew[i] = y[i] + h*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i]);
    }
}

int main() {
    double h = 0.001;
    int steps = 10000;
    int n = 2;
    

    double *t = (double *)malloc((steps+1) * sizeof(double));
    double *y1 = (double *)malloc((steps+1) * sizeof(double));
    double *y2 = (double *)malloc((steps+1) * sizeof(double));
    
    t[0] = 0.0;
    y1[0] = 1.0;
    y2[0] = 0.001;
    
    double y[2], ynew[2];
    y[0] = y1[0];
    y[1] = y2[0];
    
    for (int i = 1; i <= steps; i++) {
        rosenbrock4(t[i-1], y, h, n, ynew);
        
        t[i] = t[i-1] + h;
        y1[i] = ynew[0];
        y2[i] = ynew[1];
        
        if (isnan(y1[i]) || isnan(y2[i])) {
            printf("Обнаружены NaN значения на шаге %d. Прерывание.\n", i);
            break;
        }
        
        y[0] = y1[i];
        y[1] = y2[i];
    }
    
    FILE *fp = fopen("solution.dat", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    for (int i = 0; i <= steps; i++) {
        if (!isnan(y1[i]) && !isnan(y2[i])) {
            fprintf(fp, "%f %f %f\n", t[i], y1[i], y2[i]);
        }
    }
    
    fclose(fp);
    
    free(t);
    free(y1);
    free(y2);

    FILE *gp = fopen("plot.gp", "w");
    if (gp == NULL) {
        printf("Error opening gnuplot script file!\n");
        return 1;
    }
    
    fprintf(gp, "set terminal pngcairo enhanced font 'Arial,12'\n");
    fprintf(gp, "set output 'solution.png'\n");
    fprintf(gp, "set title 'Решение системы ОДУ методом Розенброка-Ваннера 4 порядка (α=%.1f, β=%.1f)'\n", ALPHA, BETA);
    fprintf(gp, "set xlabel 'Время t'\n");
    fprintf(gp, "set ylabel 'Значения y1, y2'\n");
    fprintf(gp, "plot 'solution.dat' using 1:2 with lines title 'y1(t)', \\\n");
    fprintf(gp, "     'solution.dat' using 1:3 with lines title 'y2(t)'\n");
    
    fclose(gp);
    

    if (system("gnuplot plot.gp") != 0) {
        printf("Ошибка при выполнении gnuplot\n");
        return 1;
    }
    
    printf("Решение успешно завершено. Результаты сохранены в solution.png\n");
    
    return 0;
}
