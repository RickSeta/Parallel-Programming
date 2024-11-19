#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>
#define N 2

void print_matrix(int n, double **m, double *o) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%8.2lf ", m[i][j]);
        }
        printf(" | %6.2lf", o[i]);
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
        m[i] = (double *) malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m[i][j] = rand() / (double) RAND_MAX * (10 + 8) - 8.0;
        }
    }

    return m;
}

double *fill_vector(int n) {
    double *v = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) v[i] = rand() / (double) RAND_MAX * (10 + 8) - 8.0;

    return v;
}

void change_rows(int n, double **m, int max_value_pos, int j) {
    for (int z = 0; z < n; z++) {
        double temp = m[max_value_pos][z];
        m[max_value_pos][z] = m[j][z];
        m[j][z] = temp;
    }
}


void order_matrix(int n, double **m) {
    int j = 0;
    int k = 0;

    for (int i = 0; i < n - 1; i++) {
        double max_value = m[i][i];
        int max_value_pos = i;
        for (j = n - 1; j > k; j--) {
            if (fabs(m[j][i]) > fabs(max_value)) {
                max_value = m[j][i];
                max_value_pos = j;
            }
        }
        if (max_value != m[i][i])
            change_rows(n, m, max_value_pos, j);
        k++;
    }
}

void gaussian_elimination(int n, double **m, double *v) {
    for (int i = 0; i < n; i++) {
        double pivot = m[i][i];
        m[i][i] = 1.0;

        for (int j = i + 1; j < n; j++) m[i][j] /= pivot;

        v[i] /= pivot;


        for (int j = i + 1; j < n; j++) {
            double const o = m[j][i];
            v[j] -= o * v[j - 1];
            for (int k = 0; k < n; k++) m[j][k] = m[j][k] - o * m[i][k];
        }

        for (int j = i - 1; j > -1; j--) {
            double const o = m[j][i];
            v[j] -= o * v[j + 1];
            for (int k = 0; k < n; k++) m[j][k] = m[j][k] - o * m[i][k];
        }
    }
}

void free_matrix(int n, double **m) {
    for (int i = 0; i < n; i++) free(m[i]);
    free(m);
}


void serial_gaussian_elimination(int n, double **m, double *o) {
    order_matrix(n, m);

    puts("Matriz inicial:\n");
    print_matrix(n, m, o);

    gaussian_elimination(n, m, o);

    puts("Matriz resultante:\n");
    print_matrix(n, m, o);

    free_matrix(n, m);

    free(o);
}

int main(int argc, char **argv) {
    int rank, procs;
    double **m, *o, start_time, end_time;

    srand(time(NULL));

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs == 1) {
        start_time = MPI_Wtime();
        m = fill_matrix(N);
        o = fill_vector(N);
        serial_gaussian_elimination(N, m, o);
        end_time = MPI_Wtime();
        printf("Tempo gasto: %lf\n", end_time - start_time);
    } else {
        int const num_lines = N / procs;

        if (rank == 0) {
            m = fill_matrix(N);
            o = fill_vector(N);
            start_time = MPI_Wtime();
            order_matrix(N, m);
        }

        for (int j = num_lines * rank; j < num_lines * rank + num_lines; j++) {
            double pivot = m[j][j];
            MPI_Bcast(&pivot, 1,MPI_DOUBLE, rank,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        if (rank == 0) {
            puts("Resultado:\n");
            print_matrix(N, m, o);
            end_time = MPI_Wtime();
            printf("Tempo gasto: %lf\n", end_time - start_time);
            free_matrix(N, m);
            free(o);
        }
    }

    MPI_Finalize();

    return 0;
}
