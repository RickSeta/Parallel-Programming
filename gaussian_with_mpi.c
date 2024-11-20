#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>
#define N 3

void print_matrix(int n, double **m) {
    int j;
    for (int i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%8.1lf ", m[i][j]);
        }
        printf(" | %6.1lf", m[i][j]);
        puts("");
    }
    puts("");
}

int check_matrix(double **m, int n) {
    return m[n][n] == 0.0;
}

double **fill_matrix(int n) {
    double **m = malloc(sizeof(double *) * n);

    for (int i = 0; i < n; i++)
        m[i] = (double *) malloc(sizeof(double) * (n + 1));

    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n + 1; j++) {
    //         m[i][j] = rand() / (double) RAND_MAX * (10 + 8) - 8.0;
    //     }
    // }

    m[0][0] = 3.0;
    m[0][1] = 2.0;
    m[0][2] = -4.0;
    m[1][0] = 2.0;
    m[1][1] = 3.0;
    m[1][2] = 3.0;
    m[2][0] = 5.0;
    m[2][1] = -3.0;
    m[2][2] = 1.0;
    m[0][3] = 3;
    m[1][3] = 15;
    m[2][3] = 14;

    return m;
}

void change_rows(int n, double **m, int max_value_pos, int j) {
    for (int z = 0; z < n + 1; z++) {
        double temp = m[max_value_pos][z];
        m[max_value_pos][z] = m[j][z];
        m[j][z] = temp;
    }
}


void order_matrix(int n, double **m, int i) {
    double max_value = m[i][i];
    int max_value_pos = i;

    for (int j = n - 1; j > -1; j--) {
        if (fabs(m[j][i]) > fabs(max_value)) {
            max_value = m[j][i];
            max_value_pos = j;
        }
    }
    if (max_value != m[i][i])
        change_rows(n, m, max_value_pos, i);
}

void gaussian_elimination(int n, double **m) {
    for (int i = 0; i < n; i++) {
        order_matrix(n, m, i);

        // print_matrix(n, m);

        double pivot = m[i][i];
        m[i][i] = 1.0;

        for (int j = i + 1; j < n + 1; j++) m[i][j] /= pivot;

        for (int j = i + 1; j < n; j++) {
            double factor = m[j][i];
            for (int k = 0; k < n + 1; k++) m[j][k] = m[j][k] - factor * m[i][k];
        }
    }

    // print_matrix(n, m);

    for (int i = n - 1; i > -1; i--) {
        for (int j = i - 1; j > -1; j--) {
            double factor = m[j][i];
            for (int k = 0; k < n + 1; k++) m[j][k] = m[j][k] - factor * m[i][k];
        }
    }
}

void free_matrix(int n, double **m) {
    for (int i = 0; i < n; i++) free(m[i]);
    free(m);
}


void serial_gaussian_elimination(int n, double **m) {
    puts("Matriz inicial:\n");
    print_matrix(n, m);

    gaussian_elimination(n, m);

    puts("Matriz resultante:\n");
    print_matrix(n, m);

    free_matrix(n, m);
}

int main(int argc, char **argv) {
    int rank, procs;
    double **m, start_time, end_time;

    srand(time(NULL));

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs == 1) {
        start_time = MPI_Wtime();
        m = fill_matrix(N);
        serial_gaussian_elimination(N, m);
        end_time = MPI_Wtime();
        printf("Tempo gasto: %lf\n", end_time - start_time);
    } else {
        int const num_lines = N / procs;

        if (rank == 0) {
            m = fill_matrix(N);
            start_time = MPI_Wtime();
        }

        for (int j = num_lines * rank; j < num_lines * rank + num_lines; j++) {
            double pivot = m[j][j];
            MPI_Bcast(&pivot, 1,MPI_DOUBLE, rank,MPI_COMM_WORLD);

        }

        if (rank == 0) {
            puts("Resultado:\n");
            print_matrix(N, m);
            end_time = MPI_Wtime();
            printf("Tempo gasto: %lf\n", end_time - start_time);
            free_matrix(N, m);
        }
    }

    MPI_Finalize();

    return 0;
}
