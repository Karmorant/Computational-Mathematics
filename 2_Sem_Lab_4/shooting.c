#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000
#define MAX_ITER 100
#define TOL 1e-6


void newton_method(double h, double* x, double* y);


double f_shooting(double x, double y) {
    return x * sqrt(y);
}

void shooting_method(double h, double* x, double* y) {
    double alpha0 = 1.0;
    double alpha1 = 2.0;
    double y1_0, y1_1;
    double u[N];
    int iter = 0;
    
    // Первый выстрел
    y[0] = 0.0;
    u[0] = alpha0;
    for (int i = 1; i < N; i++) {
        x[i] = x[i-1] + h;
        double k1 = h * u[i-1];
        double l1 = h * f_shooting(x[i-1], y[i-1]);
        
        double k2 = h * (u[i-1] + l1/2);
        double l2 = h * f_shooting(x[i-1] + h/2, y[i-1] + k1/2);
        
        double k3 = h * (u[i-1] + l2/2);
        double l3 = h * f_shooting(x[i-1] + h/2, y[i-1] + k2/2);
        
        double k4 = h * (u[i-1] + l3);
        double l4 = h * f_shooting(x[i-1] + h, y[i-1] + k3);
        
        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6;
        u[i] = u[i-1] + (l1 + 2*l2 + 2*l3 + l4)/6;
    }
    y1_0 = y[N-1];
    
    // Второй выстрел
    y[0] = 0.0;
    u[0] = alpha1;
    for (int i = 1; i < N; i++) {
        x[i] = x[i-1] + h;
        double k1 = h * u[i-1];
        double l1 = h * f_shooting(x[i-1], y[i-1]);
        
        double k2 = h * (u[i-1] + l1/2);
        double l2 = h * f_shooting(x[i-1] + h/2, y[i-1] + k1/2);
        
        double k3 = h * (u[i-1] + l2/2);
        double l3 = h * f_shooting(x[i-1] + h/2, y[i-1] + k2/2);
        
        double k4 = h * (u[i-1] + l3);
        double l4 = h * f_shooting(x[i-1] + h, y[i-1] + k3);
        
        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6;
        u[i] = u[i-1] + (l1 + 2*l2 + 2*l3 + l4)/6;
    }
    y1_1 = y[N-1];
    
    // Итерационный процесс
    while (iter < MAX_ITER && fabs(y1_1 - 2.0) > TOL) {
        double alpha_new = alpha1 - (y1_1 - 2.0) * (alpha1 - alpha0) / (y1_1 - y1_0);
        alpha0 = alpha1;
        alpha1 = alpha_new;
        y1_0 = y1_1;
        
        y[0] = 0.0;
        u[0] = alpha1;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
            double k1 = h * u[i-1];
            double l1 = h * f_shooting(x[i-1], y[i-1]);
            
            double k2 = h * (u[i-1] + l1/2);
            double l2 = h * f_shooting(x[i-1] + h/2, y[i-1] + k1/2);
            
            double k3 = h * (u[i-1] + l2/2);
            double l3 = h * f_shooting(x[i-1] + h/2, y[i-1] + k2/2);
            
            double k4 = h * (u[i-1] + l3);
            double l4 = h * f_shooting(x[i-1] + h, y[i-1] + k3);
            
            y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6;
            u[i] = u[i-1] + (l1 + 2*l2 + 2*l3 + l4)/6;
        }
        y1_1 = y[N-1];
        iter++;
    }
}

int main() {
    double h = 1.0 / (N - 1);
    double x[N], y_shooting[N], y_newton[N];
    
    // Инициализация сетки
    x[0] = 0.0;
    for (int i = 1; i < N; i++) {
        x[i] = x[i-1] + h;
    }
    

    shooting_method(h, x, y_shooting);
    

    newton_method(h, x, y_newton);
    

    printf("x\tShooting\tNewton\tDifference\n");
    for (int i = 0; i < N; i += N/10) {
        printf("%.3f\t%.6f\t%.6f\t%.6f\n", 
               x[i], y_shooting[i], y_newton[i], 
               fabs(y_shooting[i] - y_newton[i]));
    }
    
    return 0;
}