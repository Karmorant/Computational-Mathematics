#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DOMAIN_LENGTH 8.0
#define FINAL_TIME 3.0
#define WAVE_SPEED 0.5
#define COURANT_FACTOR 0.8
#define PI 3.141592653589793

// Начальное условие
void start_condition(int n, double* x, double* u) {
    for (int i = 0; i < n; i++) {
        u[i] = sin(4 * PI * x[i] / DOMAIN_LENGTH);
    }
}

// Аналитическое решение
void true_solution(int n, double* x, double t, double* u) {
    for (int i = 0; i < n; i++) {
        double pos = fmod(x[i] - WAVE_SPEED * t, DOMAIN_LENGTH);
        if (pos < 0) pos += DOMAIN_LENGTH;
        u[i] = sin(4 * PI * pos / DOMAIN_LENGTH);
    }
}

// Схема уголок
void upwind_method(int n, double* u, double speed, double dt, double dx) {
    double* u_new = malloc(n * sizeof(double));
    double factor = speed * dt / dx;
    
    for (int i = 1; i < n; i++) {
        u_new[i] = u[i] - factor * (u[i] - u[i-1]);
    }
    u_new[0] = u[0] - factor * (u[0] - u[n-1]);
    
    memcpy(u, u_new, n * sizeof(double));
    free(u_new);
}

// Схема прямоугольник
void central_diff_method(int n, double* u, double speed, double dt, double dx) {
    double* u_new = malloc(n * sizeof(double));
    double factor = speed * dt / (2 * dx);
    
    for (int i = 1; i < n-1; i++) {
        u_new[i] = u[i] - factor * (u[i+1] - u[i-1]);
    }
    u_new[0] = u[0] - factor * (u[1] - u[n-1]);
    u_new[n-1] = u[n-1] - factor * (u[0] - u[n-2]);
    
    memcpy(u, u_new, n * sizeof(double));
    free(u_new);
}

// Схема Лакс-Вендрофф
void lax_wendroff_method(int n, double* u, double speed, double dt, double dx) {
    double* u_new = malloc(n * sizeof(double));
    double sigma = speed * dt / dx;
    double sigma_sq = sigma * sigma;
    
    for (int i = 1; i < n-1; i++) {
        u_new[i] = u[i] - (sigma/2)*(u[i+1]-u[i-1]) + (sigma_sq/2)*(u[i+1]-2*u[i]+u[i-1]);
    }
    u_new[0] = u[0] - (sigma/2)*(u[1]-u[n-1]) + (sigma_sq/2)*(u[1]-2*u[0]+u[n-1]);
    u_new[n-1] = u[n-1] - (sigma/2)*(u[0]-u[n-2]) + (sigma_sq/2)*(u[0]-2*u[n-1]+u[n-2]);
    
    memcpy(u, u_new, n * sizeof(double));
    free(u_new);
}

// Сохранение результатов в файл
void save_solution(const char* filename, int n, double* x, double* u) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        perror("Failed to open file");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%f %f\n", x[i], u[i]);
    }
    fclose(fp);
}

// Построение графиков с помощью gnuplot
void plot_results(int n) {
    FILE* gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        perror("Failed to open gnuplot");
        exit(1);
    }
    fprintf(gp, "set terminal pngcairo enhanced font 'arial,10' size 1000,600\n");
    fprintf(gp, "set output 'solutions.png'\n");
    fprintf(gp, "set title 'Решения уравнения переноса'\n");
    fprintf(gp, "set xlabel 'Координата'\n");
    fprintf(gp, "set ylabel 'Амплитуда'\n");
    fprintf(gp, "plot 'true_solution.dat' w l lw 2 title 'Аналитическое', ");
    fprintf(gp, "'upwind.dat' w l lw 2 title 'Уголок', ");
    fprintf(gp, "'central.dat' w l lw 2 title 'Центральная', ");
    fprintf(gp, "'lax_wendroff.dat' w l lw 2 title 'Лакс-Вендрофф'\n");
    pclose(gp);
}

// Анализ сходимости
void convergence_analysis() {
    int node_counts[] = {40, 80, 160, 320, 640, 1280};
    int num_tests = sizeof(node_counts)/sizeof(node_counts[0]);
    
    FILE* fp = fopen("convergence.dat", "w");
    if (!fp) {
        perror("Failed to open convergence file");
        exit(1);
    }
    
    for (int test = 0; test < num_tests; test++) {
        int n = node_counts[test];
        double dx = DOMAIN_LENGTH / (n - 1);
        double dt = COURANT_FACTOR * dx / WAVE_SPEED;
        int num_steps = (int)(FINAL_TIME / dt);
        
        double* x = malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
            x[i] = i * dx;
        }
        
        double* u_true = malloc(n * sizeof(double));
        double* u_upwind = malloc(n * sizeof(double));
        double* u_central = malloc(n * sizeof(double));
        double* u_lw = malloc(n * sizeof(double));
        
        start_condition(n, x, u_upwind);
        start_condition(n, x, u_central);
        start_condition(n, x, u_lw);
        
        for (int step = 0; step < num_steps; step++) {
            upwind_method(n, u_upwind, WAVE_SPEED, dt, dx);
            central_diff_method(n, u_central, WAVE_SPEED, dt, dx);
            lax_wendroff_method(n, u_lw, WAVE_SPEED, dt, dx);
        }
        
        true_solution(n, x, FINAL_TIME, u_true);
        
        double error_upwind = 0.0, error_central = 0.0, error_lw = 0.0;
        for (int i = 0; i < n; i++) {
            error_upwind += (u_upwind[i] - u_true[i]) * (u_upwind[i] - u_true[i]);
            error_central += (u_central[i] - u_true[i]) * (u_central[i] - u_true[i]);
            error_lw += (u_lw[i] - u_true[i]) * (u_lw[i] - u_true[i]);
        }
        error_upwind = sqrt(error_upwind * dx);
        error_central = sqrt(error_central * dx);
        error_lw = sqrt(error_lw * dx);
        
        fprintf(fp, "%d %f %e %e %e\n", n, dx, error_upwind, error_central, error_lw);
        
        free(x); free(u_true); free(u_upwind); free(u_central); free(u_lw);
    }
    fclose(fp);
    
    // Построение графика сходимости
    FILE* gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        perror("Failed to open gnuplot");
        exit(1);
    }
    fprintf(gp, "set terminal pngcairo enhanced font 'arial,10' size 1000,600\n");
    fprintf(gp, "set output 'convergence.png'\n");
    fprintf(gp, "set title 'Сходимость численных методов'\n");
    fprintf(gp, "set xlabel 'Размер шага dx'\n");
    fprintf(gp, "set ylabel 'Ошибка (L2)'\n");
    fprintf(gp, "set logscale xy\n");
    fprintf(gp, "plot 'convergence.dat' u 2:3 w lp title 'Уголок', ");
    fprintf(gp, "'convergence.dat' u 2:4 w lp title 'Центральная', ");
    fprintf(gp, "'convergence.dat' u 2:5 w lp title 'Лакс-Вендрофф', ");
    fprintf(gp, "x*0.3 title 'O(dx)', x*x*3 title 'O(dx²)'\n");
    pclose(gp);
}

int main() {
    int n = 120;
    double dx = DOMAIN_LENGTH / (n - 1);
    double dt = COURANT_FACTOR * dx / WAVE_SPEED;
    int num_steps = (int)(FINAL_TIME / dt);
    
    double* x = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        x[i] = i * dx;
    }
    
    double* u_true = malloc(n * sizeof(double));
    double* u_upwind = malloc(n * sizeof(double));
    double* u_central = malloc(n * sizeof(double));
    double* u_lw = malloc(n * sizeof(double));
    
    start_condition(n, x, u_upwind);
    start_condition(n, x, u_central);
    start_condition(n, x, u_lw);
    
    for (int step = 0; step < num_steps; step++) {
        upwind_method(n, u_upwind, WAVE_SPEED, dt, dx);
        central_diff_method(n, u_central, WAVE_SPEED, dt, dx);
        lax_wendroff_method(n, u_lw, WAVE_SPEED, dt, dx);
    }
    
    true_solution(n, x, FINAL_TIME, u_true);
    
    save_solution("true_solution.dat", n, x, u_true);
    save_solution("upwind.dat", n, x, u_upwind);
    save_solution("central.dat", n, x, u_central);
    save_solution("lax_wendroff.dat", n, x, u_lw);
    
    plot_results(n);
    convergence_analysis();
    
    free(x); free(u_true); free(u_upwind); free(u_central); free(u_lw);
    return 0;
}