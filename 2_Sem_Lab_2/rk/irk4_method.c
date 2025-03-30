#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define BETA 10.0
#define ALPHA 1000.0
#define SAFE_DIV 1e-12

void higgins_system(double t, const double y[], double f[]) {

    double denom = y[1] + BETA;
    if (fabs(denom) < SAFE_DIV) {
        denom = copysign(SAFE_DIV, denom);
    }
    
    f[0] = 1.0 - y[0] * y[1];
    f[1] = ALPHA * y[1] * (y[0] - (1.0 + BETA)/denom);
}

void rk4_step(double t, double y[], double h, double y_next[]) {
    double k1[2], k2[2], k3[2], k4[2], temp[2];
    
    higgins_system(t, y, k1);
    
    temp[0] = y[0] + 0.5*h*k1[0];
    temp[1] = y[1] + 0.5*h*k1[1];
    higgins_system(t + 0.5*h, temp, k2);
    
    temp[0] = y[0] + 0.5*h*k2[0];
    temp[1] = y[1] + 0.5*h*k2[1];
    higgins_system(t + 0.5*h, temp, k3);
    
    temp[0] = y[0] + h*k3[0];
    temp[1] = y[1] + h*k3[1];
    higgins_system(t + h, temp, k4);
    
    y_next[0] = y[0] + h/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    y_next[1] = y[1] + h/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    
    if (isnan(y_next[0]) || isnan(y_next[1])) {
        y_next[0] = y[0];
        y_next[1] = y[1];
    }
}

int main() {
    double t0 = 0.0;
    double y0[2] = {1.0, 0.001};  
    double h = 0.0005;     
    int steps = 20000;             
    
    FILE *fp = fopen("higgins_rk4.dat", "w");
    if (!fp) {
        perror("Ошибка открытия файла");
        return 1;
    }
    
    double t = t0;
    double y[2] = {y0[0], y0[1]};
    
    fprintf(fp, "%.6f %.12f %.12f\n", t, y[0], y[1]);
    
    for (int i = 0; i < steps; i++) {
        double y_new[2];
        rk4_step(t, y, h, y_new);
        
        t += h;
        y[0] = y_new[0];
        y[1] = y_new[1];
        
        fprintf(fp, "%.6f %.12f %.12f\n", t, y[0], y[1]);
        
        if (isnan(y[0])) {
            printf("Ошибка: NaN на шаге %d\n", i);
            break;
        }
    }
    
    fclose(fp);
    printf("Результаты сохранены в higgins_rk4.dat\n");
    
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    fprintf(gnuplot, "set terminal pngcairo enhanced\n");
    fprintf(gnuplot, "set output 'higgins_rk4.png'\n");
    fprintf(gnuplot, "set title 'Система Хиггинса (RK4, α=%.0f)'\n", ALPHA);
    fprintf(gnuplot, "set xlabel 't'\n");
    fprintf(gnuplot, "plot 'higgins_rk4.dat' using 1:2 with lines title 'y₁(t)', ");
    fprintf(gnuplot, "'higgins_rk4.dat' using 1:3 with lines title 'y₂(t)'\n");
    pclose(gnuplot);
    
    return 0;
}