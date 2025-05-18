#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000
#define MAX_ITER 100
#define TOL 1e-6

double f_newton(double x, double y) {
    return x * sqrt(y);
}

void newton_method(double h, double* x, double* y) {
    double y_prev[N], y_new[N];
    double a[N], b[N], c[N], d[N];
    int iter = 0;
    double error;
    
    // Начальное приближение
    for (int i = 0; i < N; i++) {
        y_prev[i] = 2.0 * x[i];
    }
    
    do {
        // Построение трехдиагональной системы
        for (int i = 1; i < N-1; i++) {
            double df_dy = -0.5 * x[i] / sqrt(y_prev[i]);
            
            a[i] = 1.0 / (h*h);
            b[i] = -2.0 / (h*h) + df_dy;
            c[i] = 1.0 / (h*h);
            d[i] = f_newton(x[i], y_prev[i]) - df_dy * y_prev[i];
        }
        
        // Граничные условия
        b[0] = 1.0;
        c[0] = 0.0;
        d[0] = 0.0;
        
        a[N-1] = 0.0;
        b[N-1] = 1.0;
        d[N-1] = 2.0;
        
        // Метод прогонки
        for (int i = 1; i < N; i++) {
            double m = a[i] / b[i-1];
            b[i] = b[i] - m * c[i-1];
            d[i] = d[i] - m * d[i-1];
        }
        
        y_new[N-1] = d[N-1] / b[N-1];
        for (int i = N-2; i >= 0; i--) {
            y_new[i] = (d[i] - c[i] * y_new[i+1]) / b[i];
        }
        
        // Проверка сходимости
        error = 0.0;
        for (int i = 0; i < N; i++) {
            error += fabs(y_new[i] - y_prev[i]);
            y_prev[i] = y_new[i];
        }
        error /= N;
        
        iter++;
    } while (iter < MAX_ITER && error > TOL);
    
    // Сохраняем результат
    for (int i = 0; i < N; i++) {
        y[i] = y_new[i];
    }
}