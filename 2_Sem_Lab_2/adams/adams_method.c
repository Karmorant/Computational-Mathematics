#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define BETA 10.0
#define ALPHA 100.0
#define MIN_Y2 1e-12
#define MAX_Y2 1e+12

typedef struct {
    double y1, y2;
} State;

void system_ode(double t, const State* y, State* f) {
    double safe_y2 = (y->y2 < MIN_Y2) ? MIN_Y2 : (y->y2 > MAX_Y2) ? MAX_Y2 : y->y2;
    double denom = safe_y2 + BETA;
    
    f->y1 = 1.0 - y->y1 * safe_y2;
    f->y2 = ALPHA * safe_y2 * (y->y1 - (1.0 + BETA)/denom);
}

bool rk4_step(double t, const State* y, double h, State* y_next) {
    State k1, k2, k3, k4, temp;
    
    system_ode(t, y, &k1);
    
    temp.y1 = y->y1 + 0.5*h*k1.y1;
    temp.y2 = fmax(fmin(y->y2 + 0.5*h*k1.y2, MAX_Y2), MIN_Y2);
    system_ode(t + 0.5*h, &temp, &k2);
    
    temp.y1 = y->y1 + 0.5*h*k2.y1;
    temp.y2 = fmax(fmin(y->y2 + 0.5*h*k2.y2, MAX_Y2), MIN_Y2);
    system_ode(t + 0.5*h, &temp, &k3);
    
    temp.y1 = y->y1 + h*k3.y1;
    temp.y2 = fmax(fmin(y->y2 + h*k3.y2, MAX_Y2), MIN_Y2);
    system_ode(t + h, &temp, &k4);
    
    y_next->y1 = y->y1 + h/6.0*(k1.y1 + 2*k2.y1 + 2*k3.y1 + k4.y1);
    y_next->y2 = fmax(fmin(y->y2 + h/6.0*(k1.y2 + 2*k2.y2 + 2*k3.y2 + k4.y2), MAX_Y2), MIN_Y2);
    
    return !(isnan(y_next->y1) || isnan(y_next->y2));
}

bool adams_bashforth_4step(const State y_hist[4], const State f_hist[4], double h, State* y_next) {
    const double coeffs[4] = {55.0/24.0, -59.0/24.0, 37.0/24.0, -9.0/24.0};
    
    y_next->y1 = y_hist[3].y1 + h*(coeffs[0]*f_hist[3].y1 + coeffs[1]*f_hist[2].y1 + 
                                   coeffs[2]*f_hist[1].y1 + coeffs[3]*f_hist[0].y1);
    
    y_next->y2 = y_hist[3].y2 + h*(coeffs[0]*f_hist[3].y2 + coeffs[1]*f_hist[2].y2 + 
                                   coeffs[2]*f_hist[1].y2 + coeffs[3]*f_hist[0].y2);
    
    y_next->y2 = fmax(fmin(y_next->y2, MAX_Y2), MIN_Y2);
    
    return !(isnan(y_next->y1) || isnan(y_next->y2));
}

void solve_system(double t0, const State* y0, double h, int max_steps, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        perror("File error");
        exit(EXIT_FAILURE);
    }
    
    State y_hist[4], f_hist[4];
    double t = t0;
    
    y_hist[0] = *y0;
    system_ode(t, &y_hist[0], &f_hist[0]);
    fprintf(fp, "%.6f %.12f %.12f\n", t, y_hist[0].y1, y_hist[0].y2);
    
    // Разгонный участок
    for (int i = 1; i < 4; i++) {
        if (!rk4_step(t, &y_hist[i-1], h, &y_hist[i])) {
            fclose(fp);
            fprintf(stderr, "RK4 failed at startup step %d\n", i);
            exit(EXIT_FAILURE);
        }
        t += h;
        system_ode(t, &y_hist[i], &f_hist[i]);
        fprintf(fp, "%.6f %.12f %.12f\n", t, y_hist[i].y1, y_hist[i].y2);
    }
    
    for (int step = 4; step < max_steps; step++) {
        State y_next;
        
        if (!adams_bashforth_4step(y_hist, f_hist, h, &y_next)) {
            fprintf(stderr, "Adams-Bashforth failed at step %d (t=%.5f)\n", step, t);
            break;
        }
        
        t += h;
        system_ode(t, &y_next, &f_hist[3]);
        
        for (int i = 0; i < 3; i++) {
            y_hist[i] = y_hist[i+1];
            f_hist[i] = f_hist[i+1];
        }
        y_hist[3] = y_next;
        
        fprintf(fp, "%.6f %.12f %.12f\n", t, y_next.y1, y_next.y2);
    }
    
    fclose(fp);
}

int main() {
    State y0 = {1.0, 0.001};
    double h = 0.0075;
    int steps = 2000; 
    
    solve_system(0.0, &y0, h, steps, "adams_results.dat");
    printf("Calculation completed. Results saved to adams_results.dat\n");
    
    return EXIT_SUCCESS;
}