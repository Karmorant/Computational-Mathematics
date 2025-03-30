#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define SIGMA 10.0
#define R 28.0
#define B 8.0/3.0

void derivatives(double t, double x[], double dx[]) {
    dx[0] = -SIGMA * (x[0] - x[1]);
    dx[1] = -x[0] * x[2] + R * x[0] - x[1];
    dx[2] = x[0] * x[1] - B * x[2];
}

void rk4_step(double t, double x[], double h, double x_next[]) {
    double k1[3], k2[3], k3[3], k4[3], temp[3];
    
    derivatives(t, x, k1);
    for (int i = 0; i < 3; i++) temp[i] = x[i] + h/2 * k1[i];
    
    derivatives(t + h/2, temp, k2);
    for (int i = 0; i < 3; i++) temp[i] = x[i] + h/2 * k2[i];
    
    derivatives(t + h/2, temp, k3);
    for (int i = 0; i < 3; i++) temp[i] = x[i] + h * k3[i];
    
    derivatives(t + h, temp, k4);
    for (int i = 0; i < 3; i++)
        x_next[i] = x[i] + h/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

void bdf2(double t0, double x0[], double h, int n, FILE *out) {
    double *t = malloc((n+1)*sizeof(double));
    double (*x)[3] = malloc((n+1)*sizeof(double[3]));
    double (*dx)[3] = malloc((n+1)*sizeof(double[3]));
    
    t[0] = t0;
    for (int i = 0; i < 3; i++) x[0][i] = x0[i];
    derivatives(t[0], x[0], dx[0]);
    
    rk4_step(t[0], x[0], h, x[1]);
    t[1] = t[0] + h;
    derivatives(t[1], x[1], dx[1]);
    fprintf(out, "%f %f %f %f\n", t[1], x[1][0], x[1][1], x[1][2]);
    
    for (int step = 1; step < n; step++) {
        t[step+1] = t[step] + h;
        for (int i = 0; i < 3; i++) {
            x[step+1][i] = (4*x[step][i] - x[step-1][i] + 2*h*dx[step][i]) / 3;
        }
        derivatives(t[step+1], x[step+1], dx[step+1]);
        fprintf(out, "%f %f %f %f\n", t[step+1], x[step+1][0], x[step+1][1], x[step+1][2]);
    }
    
    free(t); free(x); free(dx);
}


void bdf3(double t0, double x0[], double h, int n, FILE *out) {
    double *t = malloc((n+2)*sizeof(double));
    double (*x)[3] = malloc((n+2)*sizeof(double[3]));
    double (*dx)[3] = malloc((n+2)*sizeof(double[3]));
    
    t[0] = t0;
    for (int i = 0; i < 3; i++) x[0][i] = x0[i];
    derivatives(t[0], x[0], dx[0]);
    
    for (int step = 0; step < 2; step++) {
        rk4_step(t[step], x[step], h, x[step+1]);
        t[step+1] = t[step] + h;
        derivatives(t[step+1], x[step+1], dx[step+1]);
        fprintf(out, "%f %f %f %f\n", t[step+1], x[step+1][0], x[step+1][1], x[step+1][2]);
    }
    
    for (int step = 2; step < n; step++) {
        t[step+1] = t[step] + h;
        for (int i = 0; i < 3; i++) {
            x[step+1][i] = (18*x[step][i] - 9*x[step-1][i] + 2*x[step-2][i] + 6*h*dx[step][i]) / 11;
        }
        derivatives(t[step+1], x[step+1], dx[step+1]);
        fprintf(out, "%f %f %f %f\n", t[step+1], x[step+1][0], x[step+1][1], x[step+1][2]);
    }
    
    free(t); free(x); free(dx);
}


void bdf4(double t0, double x0[], double h, int n, FILE *out) {
    double *t = malloc((n+3)*sizeof(double));
    double (*x)[3] = malloc((n+3)*sizeof(double[3]));
    double (*dx)[3] = malloc((n+3)*sizeof(double[3]));
    
    t[0] = t0;
    for (int i = 0; i < 3; i++) x[0][i] = x0[i];
    derivatives(t[0], x[0], dx[0]);
    
    for (int step = 0; step < 3; step++) {
        rk4_step(t[step], x[step], h, x[step+1]);
        t[step+1] = t[step] + h;
        derivatives(t[step+1], x[step+1], dx[step+1]);
        fprintf(out, "%f %f %f %f\n", t[step+1], x[step+1][0], x[step+1][1], x[step+1][2]);
    }
    
    for (int step = 3; step < n; step++) {
        t[step+1] = t[step] + h;
        for (int i = 0; i < 3; i++) {
            x[step+1][i] = (48*x[step][i] - 36*x[step-1][i] + 16*x[step-2][i] - 3*x[step-3][i] + 12*h*dx[step][i]) / 25;
        }
        derivatives(t[step+1], x[step+1], dx[step+1]);
        fprintf(out, "%f %f %f %f\n", t[step+1], x[step+1][0], x[step+1][1], x[step+1][2]);
    }
    
    free(t); free(x); free(dx);
}

int main() {
    double t0 = 0.0;
    double x0[3] = {1.0, 1.0, 1.0};
    double h = 0.01;
    int n = 5000;
    
    char filename[20];
    sprintf(filename, "bdf%d.dat", BDF_ORDER);
    FILE *out = fopen(filename, "w");
    
    switch(BDF_ORDER) {
        case 2: bdf2(t0, x0, h, n, out); break;
        case 3: bdf3(t0, x0, h, n, out); break;
        case 4: bdf4(t0, x0, h, n, out); break;
        default: fprintf(stderr, "Unsupported BDF order\n"); return 1;
    }
    
    fclose(out);
    return 0;
}