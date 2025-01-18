
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

const int N = 4;

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

void print_array(float *array, int size, const char *label)
{
    printf("%s: [", label);
    for (int i = 0; i < size; i++)
    {
        printf("%f", array[i]);
        if (i < size - 1)
        {
            printf(", ");
        }
    }
    printf("]\n");
}

void gauss_elimination(float **A, float *x, int rank, int num_procs)
{
    // 0 > 1 ;  1 > 2;  2 > 3; 3 > 4;

    MPI_Status status;
    MPI_Request request;

    int norm, row, col, max_row;
    float multiplier;
    int i, j, k = 0;

    for (norm = 0; norm < N - 1; norm++)
    {
        MPI_Bcast(A[norm], N + 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        print_matrix(A, "matrix");
        for (row = norm + 1 + rank; row < N; row += num_procs)
        {

            // printf("\nRank %d dentro do interleaving com a row %d e norm %d e k %d", rank, row, norm, k);
            k++;
            if (fabs(A[norm][norm]) < 1e-6)
            {
                fprintf(stderr, "a matriz é singular ou quase singular.\n");
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
            multiplier = A[row][norm] / A[norm][norm];

            // print_array(A[norm], N + 1, "\n\nArray norm antes");
            // print_array(A[row], N + 1, "Array row antes");

            for (col = norm; col < N + 1; col++)
            {
                A[row][col] -= A[norm][col] * multiplier;
            }
            // printf("\nRank %d", rank);
            // print_array(A[norm], N + 1, "\nArray norm depois");
            // print_array(A[row], N + 1, "Array row depois");
        }
    }
}

void generate_matrix(float **A)
{
    int i, j, e = 0;
    srand(time(NULL));
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N + 1; j++)
        {
            A[i][j] = (rand() / 50000 + 1);
        }
    }
}

void manual_matrix(float **A)
{
    float matrix[4][5] = {
        {2, 1, -1, 3, 10},  // Row 1
        {4, -2, 1, 1, 8},   // Row 2
        {-1, 5, 3, -2, -5}, // Row 3
        {3, -3, -2, 4, 6}   // Row 4
    };

    // Copy data to the provided matrix
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            A[i][j] = matrix[i][j];
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

    print_matrix(A, "matrix do back");
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

    A = (float **)calloc(N, sizeof(float *));
    for (int q = 0; q < N; q++)
        A[q] = (float *)calloc(N + 1, sizeof(float *));

    b = (float *)malloc(sizeof(float) * N);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (rank == 0)
    {
        // generate_matrix(A);
        //  fill_diagonal_matrix(N, A);
        manual_matrix(A);
        // print_matrix(A, "matrix");

        printf("\nStart  MPI...Process=%d\n", rank);
        start_time = MPI_Wtime();
    }

    gauss_elimination(A, b, rank, num_procs);

    if (rank == 0)
    {
        back_substitution(A, b);

        end_time = MPI_Wtime();
        print_output(b);
        printf("\nelapsed time = %f\n", end_time - start_time);
    }

    for (int i = 0; i < N; i++)
        free(A[i]);
    free(b);

    MPI_Finalize();
    return 0;
}