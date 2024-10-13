#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100

void printMatrix(double mat[N][N + 1])
{
        for (int i = 0; i < N; i++)
        {
                for (int j = 0; j < N + 1; j++)
                {
                        printf("%10.4lf ", mat[i][j]);
                }
                printf("\n");
        }
        printf("\n");
}

// Функция решения СЛАУ методом Гаусса с выбором главного элемента
void gauss(double mat[N][N + 1])
{
        for (int i = 0; i < N; i++)
        {
                // Поиск главного элемента по столбцу
                int maxRow = i;
                for (int k = i + 1; k < N; k++)
                {
                        if (fabs(mat[k][i]) > fabs(mat[maxRow][i]))
                        {
                                maxRow = k;
                        }
                }

                // Поменять местами текущую строку и строку с главным элементом
                for (int k = 0; k < N + 1; k++)
                {
                        double temp = mat[maxRow][k];
                        mat[maxRow][k] = mat[i][k];
                        mat[i][k] = temp;
                }

                // Прямой ход метода Гаусса
                for (int k = i + 1; k < N; k++)
                {
                        double factor = mat[k][i] / mat[i][i];
                        for (int j = i; j < N + 1; j++)
                        {
                                mat[k][j] -= factor * mat[i][j];
                        }
                }
        }


        // Обратный ход
        double x[N];
        for (int i = N - 1; i >= 0; i--)
        {
                x[i] = mat[i][N] / mat[i][i];
                for (int k = i - 1; k >= 0; k--)
                {
                        mat[k][N] -= mat[k][i] * x[i];
                }
        }

        printf("Решения:\n");
        for (int i = 0; i < N; i++)
        {
                printf("x[%d] = %10.10lf\n", i + 1, x[i]);
        }
}

int main()
{
        double mat[N][N + 1];

        for (int i = 1; i <= N; i++)
        {
                for (int g = 1; g <= N + 1; g++)
                {
                        if (g == (N + 1))
                        {
                                mat[i - 1][g - 1] = i + 2;
                        }
                        else if (g > i + 1)
                        {
                                mat[i - 1][g - 1] = 0;
                        } 
                        else if (i == g)
                        {
                                mat[i - 1][g - 1] = 10;
                        }
                        else 
                        {
                                mat[i - 1][g - 1] = (double) g/(i + g);
                        }
                }
        }

        gauss(mat);

        return 0;
}
