// Функция для вычисления обратной матрицы
int inverseMatrix(double mat[N][N], double inv[N][N]) {
    double augMat[N][2 * N];
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            augMat[i][j] = mat[i][j];
            augMat[i][j + N] = (i == j) ? 1 : 0; 
        }
    }


    for (int i = 0; i < N; i++) {
        if (augMat[i][i] == 0) {
            printf("Обратная матрица не существует.\n");
            return 0;
        }

        double pivot = augMat[i][i];
        

        for (int j = 0; j < 2 * N; j++) {
            augMat[i][j] /= pivot;
        }


        for (int k = 0; k < N; k++) {
            if (k != i) {
                double factor = augMat[k][i];
                for (int j = 0; j < 2 * N; j++) {
                    augMat[k][j] -= factor * augMat[i][j];
                }
            }
        }
    }


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            inv[i][j] = augMat[i][j + N];
        }
    }

    return 1;
}

void addMatrices(double mat1[N][N], double mat2[N][N], double result[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
}

void multiplyMatrices(double mat1[N][N], double mat2[N][N], double result[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i][j] = 0;
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

void multiplyMatrixByScalar(double mat[N][N], double scalar, double result[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i][j] = mat[i][j] * scalar; 
        }
    }
}

void multiplyMatrixByVector(double mat[N][N], double vec[N], double result[N]) {
    for (int i = 0; i < N; i++) {
        result[i] = 0;
    }


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i] += mat[i][j] * vec[j]; 
        }
    }
}

void addVectors(double vec1[N], double vec2[N], double result[N]) {
    for (int i = 0; i < N; i++) {
        result[i] = vec1[i] + vec2[i]; 
    }
}