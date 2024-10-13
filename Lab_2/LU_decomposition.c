#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100

// Функция для LU-разложения матрицы A на L и U
void luDecomposition(double A[N][N], double L[N][N], double U[N][N])
{
        for (int i = 0; i < N; i++)
        {
                for (int k = i; k < N; k++)
                {
                        double sum = 0;
                        for (int j = 0; j < i; j++)
                                sum += (L[i][j] * U[j][k]);
                        U[i][k] = A[i][k] - sum;
                }

                for (int k = i; k < N; k++)
                {
                        if (i == k)
                                L[i][i] = 1;
                        else
                        {
                                double sum = 0;
                                for (int j = 0; j < i; j++)
                                        sum += (L[k][j] * U[j][i]);
                                L[k][i] = (A[k][i] - sum) / U[i][i];
                        }
                }
        }
}

// Функция для решения уравнения L*y = b методом прямой подстановки
void forwardSubstitution(double L[N][N], double b[N], double y[N])
{
        for (int i = 0; i < N; i++)
        {
                y[i] = b[i];
                for (int j = 0; j < i; j++)
                {
                        y[i] -= L[i][j] * y[j];
                }
        }
}

// Функция для решения уравнения U*x = y методом обратной подстановки
void backwardSubstitution(double U[N][N], double y[N], double x[N])
{
        for (int i = N - 1; i >= 0; i--)
        {
                x[i] = y[i];
                for (int j = i + 1; j < N; j++)
                {
                        x[i] -= U[i][j] * x[j];
                }
                x[i] /= U[i][i];
        }
}

// Основная функция для решения СЛАУ методом LU-разложения
void solveLU(double A[N][N], double b[N], double x[N])
{
        double L[N][N] = {0}, U[N][N] = {0}; 
        double y[N] = {0};

        luDecomposition(A, L, U);

        forwardSubstitution(L, b, y);

        backwardSubstitution(U, y, x);
}

int main()
{
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

        double x[N];

        solveLU(A, b, x);

        printf("Решение системы: \n");
        for (int i = 0; i < N; i++)
        {
                printf("x[%d] = %.10lf\n", i, x[i]);
        }

        return 0;
}