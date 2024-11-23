#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define N 4 // Size of the matrix (N x N)
#define MAX_VAL 10 // Maximum value for the random elements in the matrix

void print_matrix(double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void gaussian_elimination(double matrix[N][N], int rank, int size) {
    int i, j, k;
    double factor;

    for (i = 0; i < N; i++) {
        // Broadcast the pivot row to all processes
        MPI_Bcast(matrix[i], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Perform the row elimination for the rows below the pivot
        for (j = i + 1 + rank; j < N; j += size) {
            if (matrix[i][i] != 0) {
                factor = matrix[j][i] / matrix[i][i];
                for (k = i; k < N; k++) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
            }
        }

        // Synchronize after each step
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void back_substitution(double matrix[N][N]) {
    int i, j;
    double x[N];

    // Perform back substitution to find the solution
    for (i = N - 1; i >= 0; i--) {
        x[i] = matrix[i][N];
        for (j = i + 1; j < N; j++) {
            x[i] -= matrix[i][j] * x[j];
        }
        x[i] /= matrix[i][i];
    }

    printf("Solution: \n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %lf\n", i, x[i]);
    }
}

void generate_random_matrix(double matrix[N][N]) {
    srand(time(0));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = rand() % MAX_VAL + 1; // Random values between 1 and MAX_VAL
        }
    }
}

int main(int argc, char** argv) {
    int rank, size;
    double matrix[N][N];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        generate_random_matrix(matrix);
        printf("Generated Matrix:\n");
        print_matrix(matrix);
    }

    // Perform Gaussian elimination
    gaussian_elimination(matrix, rank, size);

    // Root process prints the final result after back substitution
    if (rank == 0) {
        back_substitution(matrix);
    }

    MPI_Finalize();
    return 0;
}
