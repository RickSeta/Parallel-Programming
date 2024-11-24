#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>
#define N 3

void print_matrix(int n, int o, double m[n][o]) {
    int j;

    // if (n <= 10) {
        for (int i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%8.3lf ", m[i][j]);
            }
            printf(" | %6.3lf", m[i][j]);
            puts("");
        }
        puts("");
    // } else {
    //     for (int i = 0; i < n; i++)
    //         printf("%8.3lf\n", m[i][n]);
    //     puts("");
    // }
}

void fill_matrix(int n, int o, double m[n][o]) {
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < o; j++) {
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
    m[0][3] = 3.0;
    m[1][3] = 15.0;
    m[2][3] = 14.0;

    // m[0][0] = 8.436; m[0][1] = 2.883; m[0][2] = -0.524; m[0][3] = -6.264;
    // m[1][0] = 3.608; m[1][1] = -6.019; m[1][2] = 4.509; m[1][3] = -7.003;
    // m[2][0] = 0.139; m[2][1] = 4.916; m[2][2] = 4.798; m[2][3] = 7.387;
}

void change_rows(int n, int o, double m[n][o], int max_value_pos, int j) {
    for (int z = 0; z < o; z++) {
        double temp = m[max_value_pos][z];
        m[max_value_pos][z] = m[j][z];
        m[j][z] = temp;
    }
}


void order_matrix(int n, int o, double m[n][o], int i) {
    double max_value = m[i][i];
    int max_value_pos = i;

    for (int j = n - 1; j > -1; j--) {
        if (fabs(m[j][i]) > fabs(max_value)) {
            max_value = m[j][i];
            max_value_pos = j;
        }
    }
    if (max_value != m[i][i])
        change_rows(n, n + 1, m, max_value_pos, i);
}


void gaussian_elimination(int n, int o, double m[n][o]) {
    for (int i = 0; i < n; i++) {
        order_matrix(n, o, m, i);

        double pivot = m[i][i];
        m[i][i] = 1.0;

        for (int j = i + 1; j < n + 1; j++) m[i][j] /= pivot;

        for (int j = i + 1; j < n; j++) {
            double factor = m[j][i];
            for (int k = 0; k < o; k++) m[j][k] = m[j][k] - factor * m[i][k];
        }
    }

    for (int i = n - 1; i > -1; i--) {
        for (int j = i - 1; j > -1; j--) {
            double factor = m[j][i];
            for (int k = 0; k < o; k++) m[j][k] = m[j][k] - factor * m[i][k];
        }
    }
}

void serial_gaussian_elimination(int n, int o, double m[n][o]) {
    puts("Valores Iniciais:\n");
    print_matrix(n, o, m);

    gaussian_elimination(n, o, m);

    puts("Valores Resultantes:\n");
    print_matrix(n, o, m);
}

void fill_pivot_line(int n, int o, double pivot_line[o], double m[n][o], int i) {
    for (int z = 0; z < o; z++) {
        pivot_line[z] = m[i][z];
    }
}

void back_elimination(int n, int o, double m[n][o]) {
    for (int i = n - 1; i > -1; i--) {
        for (int j = i - 1; j > -1; j--) {
            double factor = m[j][i];
            for (int k = 0; k < o; k++) m[j][k] = m[j][k] - factor * m[i][k];
        }
    }
}

int main(int argc, char **argv) {
    int rank, procs;
    double m[N][N + 1], start_time, end_time;

    srand(time(NULL));

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);


    // algoritmo sequencial
    if (procs == 1) {
        start_time = MPI_Wtime();
        fill_matrix(N, N + 1, m);
        serial_gaussian_elimination(N, N + 1, m);
        end_time = MPI_Wtime();
        printf("Tempo gasto: %lf\n", end_time - start_time);
    }
    // algoritmo paralelo
    else {
        if (N % procs != 0 ) {
            puts("O número de linhas deve ser igual ou maior ao de processos");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // linhas por processos
        int num_lines = N / procs;

        double local_matrix[num_lines][N + 1];

        double pivot_line[N + 1];


        if (rank == 0) {
            start_time = MPI_Wtime();
            fill_matrix(N,N + 1, m);
            print_matrix(N,N + 1, m);
        }

        for (int i = 0; i < N; i++) {
            if (rank == 0) {
                order_matrix(N, N + 1, m, i);
                fill_pivot_line(N, N + 1, pivot_line, m, i);
                // double pivot = pivot_line[i];
                // for (int j = 0; j < N + 1; j++) {
                //     pivot_line[j] /= pivot;
                //     printf("%lf ", pivot_line[j]);
                // }
                // puts("");
            }

            MPI_Scatter(*m,
                        (N + 1) * num_lines,
                        MPI_DOUBLE,
                        *local_matrix,
                        (N + 1) * num_lines,
                        MPI_DOUBLE, 0,
                        MPI_COMM_WORLD);

            MPI_Bcast(pivot_line, N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (i / num_lines == rank) {
                for (int j = i; j < N + 1; j++) local_matrix[i % num_lines][j] /= pivot_line[i];

                for (int j = i % num_lines + 1; j < num_lines; j++) {
                    double factor = local_matrix[j][i];
                    for (int k = i; k < N + 1; k++) local_matrix[j][k] -= factor * pivot_line[k] / pivot_line[i];
                }
            }
            else {
                for (int j = 0; j < num_lines; j++) {
                    double factor = local_matrix[j][i];
                    for (int k = i; k < N + 1; k++) local_matrix[j][k] -= factor * pivot_line[k] / pivot_line[i];
                }
            }


            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Gather(*local_matrix,
                       (N + 1) * num_lines,
                       MPI_DOUBLE,
                       *m,
                       (N + 1) * num_lines,
                       MPI_DOUBLE, 0,
                       MPI_COMM_WORLD);
        }

        MPI_Gather(*local_matrix,
                   (N + 1) * num_lines,
                   MPI_DOUBLE,
                   *m,
                   (N + 1) * num_lines,
                   MPI_DOUBLE, 0,
                   MPI_COMM_WORLD);

        // resultado final
        if (rank == 0) {
            back_elimination(N, N+1, m);
            puts("Resultado:\n");
            print_matrix(N,N + 1, m);
            end_time = MPI_Wtime();
            printf("Tempo gasto: %lf\n", end_time - start_time);
        }
    }

    MPI_Finalize();

    return 0;
}

// // rank de onde vem o pivot
// double max_pivot_location;
//
// // criacao de sub matriz
// double local_matrix[num_lines][N + 1];
//
// //linhas pivot
// double local_pivot_line[N + 2], global_pivot_line[N + 2];
// double all_pivot_lines[procs][N + 2];


// // distribui as sub matrizes
// MPI_Scatter(*m,
//             (N + 1) * num_lines,
//             MPI_DOUBLE,
//             *local_matrix,
//             (N + 1) * num_lines,
//             MPI_DOUBLE, 0,
//             MPI_COMM_WORLD);
//
//
// for (int i = rank * num_lines; i < rank * num_lines + num_lines; i++) {
//     // busca linha com maior valor
//     order_matrix(num_lines, N + 1, local_matrix, i);
//
//     //obtencao do pivot local
//     local_pivot(N, N + 1, rank, local_pivot_line, local_matrix);
//
//     MPI_Gather(local_pivot_line, procs, MPI_DOUBLE,
//                *all_pivot_lines, procs, MPI_DOUBLE,
//                rank, MPI_COMM_WORLD);
//
// }
//
// puts("E");
//
//
// // recebe os resultados de cada sub matriz
// MPI_Gather(*local_matrix,
//            (N + 1) * num_lines,
//            MPI_DOUBLE,
//            *m,
//            (N + 1) * num_lines,
//            MPI_DOUBLE, 0,
//            MPI_COMM_WORLD);


// void local_pivot(int n, int o, int r, double local_pivot[n], double local_matrix[n][o]) {
//     int j;
//     for (j = 0; j < o; j++) {
//         local_pivot[j] = local_matrix[0][j];
//     }
//     local_pivot[j] = r;
// }
//
// void find_global_pivot(int n, int o, int p, double local_pivot[n], double local_matrix[n][o]) {
// }
