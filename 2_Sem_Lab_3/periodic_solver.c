#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 200
#define PI 3.14159265358979323846

double P_squared(double x) {
    return 10.0 + sin(2 * PI * x);
}

double f_func(double x) {
    return cos(2 * PI * x);
}

void solve_periodic_tridiagonal(int n, double *a, double *b, double *c, double *d, double *x) {
    double *v = (double *)malloc(n * sizeof(double));
    double *u = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));
    
    v[0] = a[0];
    u[0] = c[n-1];
    y[0] = d[0];
    
    for (int i = 1; i < n-1; i++) {
        double m = a[i] / b[i-1];
        b[i] = b[i] - m * c[i-1];
        d[i] = d[i] - m * d[i-1];
        v[i] = -m * v[i-1];
        u[i] = (i == n-2) ? -m * u[i-1] : 0.0;
        y[i] = d[i];
    }
    
    double A = a[n-1];
    double B = b[n-1];
    double C = c[n-2];
    double D = d[n-1];
    
    for (int i = 0; i < n-2; i++) {
        A += 0;
        B += v[i] * (i == n-2 ? C : 0);
        D += y[i] * (i == n-2 ? C : 0);
    }
    
    x[n-1] = (D - A * x[0] - B * x[n-2]) / (b[n-1] + A * u[0] + B * v[n-2]);
    
    for (int i = n-2; i >= 0; i--) {
        x[i] = (d[i] - c[i] * x[i+1] - v[i] * x[n-1]) / b[i];
    }
    
    free(v);
    free(u);
    free(y);
}

int main() {
    double h = 1.0 / N;
    double x[N], a[N], b[N], c[N], d[N], y[N];
    
    for (int i = 0; i < N; i++) {
        x[i] = i * h;
        
        a[i] = 1.0 / (h * h); 
        b[i] = -2.0 / (h * h) - P_squared(x[i]);
        c[i] = 1.0 / (h * h);
        d[i] = f_func(x[i]);
    }
    
    a[0] = 0.0;
    c[N-1] = 0.0;
    
    solve_periodic_tridiagonal(N, a, b, c, d, y);
    
    FILE *fp = fopen("periodic_solution.dat", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        fprintf(fp, "%f %f\n", x[i], y[i]);
    }
    fprintf(fp, "%f %f\n", x[0], y[0]);
    
    fclose(fp);
    
    FILE *gp = fopen("plot_periodic.gp", "w");
    if (gp == NULL) {
        printf("Error opening gnuplot script file!\n");
        return 1;
    }
    
    fprintf(gp, "set terminal pngcairo enhanced font 'Arial,12'\n");
    fprintf(gp, "set output 'periodic_solution.png'\n");
    fprintf(gp, "set title 'Периодическое решение уравнения y'' - (10 + sin(2πx))y = cos(2πx)'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y(x)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot 'periodic_solution.dat' with lines title 'Решение'\n");
    
    fclose(gp);
    
    system("gnuplot plot_periodic.gp");
    
    printf("Решение успешно завершено. Результаты сохранены в periodic_solution.png\n");
    
    return 0;
}