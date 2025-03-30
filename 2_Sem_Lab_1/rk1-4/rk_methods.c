#include <stdio.h>
#include <math.h>

#define SIGMA 10.0
#define R 28.0
#define B 8.0/3.0

void derivatives(double t, double x[], double dx[]) {
    dx[0] = -SIGMA * (x[0] - x[1]);
    dx[1] = -x[0] * x[2] + R * x[0] - x[1];
    dx[2] = x[0] * x[1] - B * x[2];
}

void rk1(double t0, double x0[], double h, int n, FILE *out) {
    double x[3], dx[3];
    for (int i = 0; i < 3; i++) x[i] = x0[i];
    
    for (int step = 0; step < n; step++) {
        double t = t0 + step * h;
        derivatives(t, x, dx);
        for (int i = 0; i < 3; i++) x[i] += h * dx[i];
        fprintf(out, "%f %f %f %f\n", t+h, x[0], x[1], x[2]);
    }
}

void rk2(double t0, double x0[], double h, int n, FILE *out) {
    double x[3], k1[3], k2[3], temp[3];
    for (int i = 0; i < 3; i++) x[i] = x0[i];
    
    for (int step = 0; step < n; step++) {
        double t = t0 + step * h;
        
        derivatives(t, x, k1);
        for (int i = 0; i < 3; i++) temp[i] = x[i] + h/2 * k1[i];
        
        derivatives(t + h/2, temp, k2);
        for (int i = 0; i < 3; i++) x[i] += h * k2[i];
        
        fprintf(out, "%f %f %f %f\n", t+h, x[0], x[1], x[2]);
    }
}

void rk3(double t0, double x0[], double h, int n, FILE *out) {
    double x[3], k1[3], k2[3], k3[3], temp[3];
    for (int i = 0; i < 3; i++) x[i] = x0[i];
    
    for (int step = 0; step < n; step++) {
        double t = t0 + step * h;
        
        derivatives(t, x, k1);
        for (int i = 0; i < 3; i++) temp[i] = x[i] + h/2 * k1[i];
        
        derivatives(t + h/2, temp, k2);
        for (int i = 0; i < 3; i++) temp[i] = x[i] + h * (-k1[i] + 2*k2[i]);
        
        derivatives(t + h, temp, k3);
        for (int i = 0; i < 3; i++) 
            x[i] += h/6 * (k1[i] + 4*k2[i] + k3[i]);
        
        fprintf(out, "%f %f %f %f\n", t+h, x[0], x[1], x[2]);
    }
}

void rk4(double t0, double x0[], double h, int n, FILE *out) {
    double x[3], k1[3], k2[3], k3[3], k4[3], temp[3];
    for (int i = 0; i < 3; i++) x[i] = x0[i];
    
    for (int step = 0; step < n; step++) {
        double t = t0 + step * h;
        
        derivatives(t, x, k1);
        for (int i = 0; i < 3; i++) temp[i] = x[i] + h/2 * k1[i];
        
        derivatives(t + h/2, temp, k2);
        for (int i = 0; i < 3; i++) temp[i] = x[i] + h/2 * k2[i];
        
        derivatives(t + h/2, temp, k3);
        for (int i = 0; i < 3; i++) temp[i] = x[i] + h * k3[i];
        
        derivatives(t + h, temp, k4);
        for (int i = 0; i < 3; i++)
            x[i] += h/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        
        fprintf(out, "%f %f %f %f\n", t+h, x[0], x[1], x[2]);
    }
}

int main() {
    double t0 = 0.0;
    double x0[3] = {1.0, 1.0, 1.0};
    double h = 0.01;
    int n = 5000;
    
    char filename[20];
    sprintf(filename, "rk%d.dat", RK_ORDER);
    FILE *out = fopen(filename, "w");
    
    switch(RK_ORDER) {
        case 1: rk1(t0, x0, h, n, out); break;
        case 2: rk2(t0, x0, h, n, out); break;
        case 3: rk3(t0, x0, h, n, out); break;
        case 4: rk4(t0, x0, h, n, out); break;
        default: fprintf(stderr, "Unsupported RK order\n"); return 1;
    }
    
    fclose(out);
    return 0;
}
