#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define N 100

void sorMethod(double A[N][N], double b[N], double x[N], double epsilon, double omega, int *iter, FILE *fp)
{
    double x_new[N];
    double resid[N];
    bool converged = false;

    *iter = 0;

    while (!converged)
    {
        for (int i = 0; i < N; i++)
        {
            double sum = 0;
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x_new[j];
                }
            }
            x_new[i] = (1 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sum);
        }

        // Проверка сходимости
        double error = 0;
        double norm = 0;
        for (int i = 0; i < N; i++)
        {
            error += fabs(x_new[i] - x[i]);
        }

        for (int i = 0; i < N; i++)
        {
            resid[i] = b[i];
            for (int j = 0; j < N; j++)
            {
                resid[i] -= A[i][j] * x_new[j];
            }
            norm += fabs(resid[i]);
        }
        fprintf(fp, "%d %.10lf\n", *iter, norm);

        if (error < epsilon)
        {
            converged = true;
        }

        for (int i = 0; i < N; i++)
        {
            x[i] = x_new[i];
        }

        (*iter)++;
    }
}

int main()
{
    FILE *fp = fopen("SOR_data.txt", "w");

    double A[N][N];
    for (int i = 1; i <= N; i++)
    {
        for (int g = 1; g <= N; g++)
        {
            if (g > i + 1)
            {
                A[i - 1][g - 1] = 0;
            }
            else if (i == g)
            {
                A[i - 1][g - 1] = 10;
            }
            else
            {
                A[i - 1][g - 1] = (double)g / (i + g);
            }
        }
    }

    double b[N];
    for (int i = 0; i < N; i++)
    {
        b[i] = i + 3;
    }

    double x[N] = {0};
    double epsilon = 1e-9;
    int iter = 0;
    double omega = 1.25;

    sorMethod(A, b, x, epsilon, omega, &iter, fp);

    printf("Решение системы уравнений:\n");
    for (int i = 0; i < N; i++)
    {
        printf("x[%d] = %.6lf\n", i, x[i]);
    }

    printf("Количество итераций: %d\n", iter);

    fclose(fp);

        FILE *gnuplotPipe = popen("./plot_2.pgi", "w");
        fclose(gnuplotPipe);

    return 0;
}
