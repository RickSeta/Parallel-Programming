#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 4  // Dimension of the matrix (N x N)
#define ROOT 0  // Root process

// Function to print the matrix
void print_matrix(int *matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%d ", matrix[i * cols + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    int *matrix = NULL;
    int *sub_matrix = NULL;
    int block_size = N * N / 4;  // For 4 processes, each block will be N/2 x N/2
    int block_rows = N / 2;
    int block_cols = N / 2;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Only the root process creates the full matrix
    if (rank == ROOT) {
        // Allocate and initialize the matrix
        matrix = (int *)malloc(N * N * sizeof(int));
        int value = 1;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                matrix[i * N + j] = value++;
            }
        }
        printf("Full Matrix on Root Process:\n");
        print_matrix(matrix, N, N);
    }

    // Allocate memory for each process to store a sub-matrix (block)
    sub_matrix = (int *)malloc(block_rows * block_cols * sizeof(int));

    // Distribute the blocks to each process
    MPI_Scatter(matrix, block_size, MPI_INT, sub_matrix, block_size, MPI_INT, ROOT, MPI_COMM_WORLD);

    // Each process computes on its block (sub-matrix)
    printf("Process %d received block:\n", rank);
    print_matrix(sub_matrix, block_rows, block_cols);

    // After computation, gather the blocks back to the root process
    MPI_Gather(sub_matrix, block_size, MPI_INT, matrix, block_size, MPI_INT, ROOT, MPI_COMM_WORLD);

    // The root process prints the final matrix after gathering the results
    if (rank == ROOT) {
        printf("Matrix after gathering blocks:\n");
        print_matrix(matrix, N, N);
        free(matrix);
    }

    free(sub_matrix);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
