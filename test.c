
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

const int N = 4;

void gauss_elimination(float **A, float *x, int rank, int num_procs)
{
    // 0 > 1 ;  1 > 2;  2 > 3; 3 > 4;

    MPI_Status status;
    MPI_Request request;

    int norm, row, col;
    float multiplier;
    int i, j;

    for (norm = 0; norm < N - 1; norm++)
    {
        MPI_Bcast(A[norm], N + 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

        for (row = norm + 1 + rank; row < N; row += num_procs)
        {
            printf("meu rank é e %i a linha que peguei é %i \n", rank, row-1);
            multiplier = A[row][norm] / A[norm][norm];
            for (col = norm; col < N + 1; col++)
            {
                A[row][col] -= A[norm][col] * multiplier;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void generate_matrix(float **A)
{
    int i, j, e = 0;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N + 1; j++)
        {
            A[i][j] = (rand() / 50000);
        }
    }
}

void fill_diagonal_matrix(int N, float **matrix)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= N; j++)
        {
            if (j == N)
            {
                matrix[i][j] = 4;
            }
            else
            {
                matrix[i][j] = (i == j) * (i + 1);
            }
        }
    }
}

void print_matrix(float **A, char *msg)
{
    int i, j;
    printf("%s\n", msg);
    for (i = 0; i < N; i++)
    {
        printf("%d\t", i);
        for (j = 0; j < N + 1; j++)
        {
            printf("%lf\t", A[i][j]);
        }
        printf("\n");
    }
}

void print_output(float *x)
{
    int i;
    for (i = 0; i < N; i++)
    {
        printf("X[%d]-->%f\n", i, x[i]);
    }
}

void back_substitution(float **A, float *x)
{
    double sum;
    int i, j;
    x[N - 1] = A[N - 1][N] / A[N - 1][N - 1];

    for (i = N - 2; i >= 0; i--)
    {
        sum = 0;
        for (j = i + 1; j <= N - 1; j++)
        {
            sum = sum + A[i][j] * x[j];
        }
        x[i] = (A[i][N] - sum) / A[i][i];
    }
}

int main(int argc, char **argv)
{
    double start_time, end_time;
    float **A, *b;
    int num_procs;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    A = (float **)calloc(N, sizeof(float *));
    for (int q = 0; q < N; q++)
        A[q] = (float *)calloc(N + 1, sizeof(float *));

    b = (float *)malloc(sizeof(float) * N);

    if (rank == 0)
    {
        // generate_matrix(A);
        fill_diagonal_matrix(N, A);

        printf("\nStart  MPI...Process=%d\n", rank);
        start_time = MPI_Wtime();
    }

    gauss_elimination(A, b, rank, num_procs);

    if (rank == 0)
    {
        back_substitution(A, b);

        print_output(b);

        end_time = MPI_Wtime();
        printf("\nelapsed time = %f\n", end_time - start_time);
    }

    for (int i = 0; i < N; i++)
        free(A[i]);
    free(b);

    MPI_Finalize();
    return 0;
}