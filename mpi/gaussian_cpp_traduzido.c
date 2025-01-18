#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>
#define N 100

void print_matrix(int n, int o, double m[n][o]) {
    int j;

    if (n <= 10) {
        for (int i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%8.3lf ", m[i][j]);
            }
            printf(" | %6.3lf", m[i][j]);
            puts("");
        }
        puts("");
    } else {
        for (int i = 0; i < n; i++)
            printf("%8.3lf\n", m[i][n]);
        puts("");
    }
}

void fill_matrix(int n, int o, double m[n][o]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < o; j++) {
            m[i][j] = rand() / (double) RAND_MAX * (10 + 8) - 8.0;
        }
    }

    // m[0][0] = 3.0;
    // m[0][1] = 2.0;
    // m[0][2] = -4.0;
    // m[1][0] = 2.0;
    // m[1][1] = 3.0;
    // m[1][2] = 3.0;
    // m[2][0] = 5.0;
    // m[2][1] = -3.0;
    // m[2][2] = 1.0;
    // m[0][3] = 3;
    // m[1][3] = 15;
    // m[2][3] = 14;

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

void back_elimination(int n, int o, double m[n][o]) {
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
        if (N % procs != 0) {
            puts("O número de linhas deve ser igual ou maior ao de processos");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // linhas por processos
        int num_lines = N / procs;
        int start_line = rank * num_lines;
        int end_line = start_line + num_lines;

        // rank de onde vem o pivot
        double max_pivot_location;

        // criacao de sub matriz
        double local_matrix[num_lines][N + 1];

        //linhas pivot
        double pivot_line[N + 1];

        MPI_Request req;

        if (rank == 0) {
            fill_matrix(N,N + 1, m);
            start_time = MPI_Wtime();
        }

        // distribui as sub matrizes
        MPI_Scatter(*m,
                    (N + 1) * num_lines,
                    MPI_DOUBLE,
                    *local_matrix,
                    (N + 1) * num_lines,
                    MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);


        // Performance gaussian elimination
        for (int lines = 0; lines < end_line; lines++) {
            // See if this process is responsible for the pivot calculation
            int check_rank = lines / num_lines;

            // If the row is mapped to this rank...
            if (rank == check_rank) {
                // Calculate the row in the local matrix
                int local_line = lines % num_lines;

                // Get the value of the pivot
                double pivot = local_matrix[local_line][lines];

                // Divide the rest of the row by the pivot
                for (int col = lines; col < N + 1; col++) {
                    local_matrix[local_line][col] /= pivot;
                }

                // Send the pivot row to the other processes
                for (int i = check_rank + 1; i < procs; i++) {
                    MPI_Isend(&local_matrix[local_line][0], N, MPI_DOUBLE, i, 0,
                              MPI_COMM_WORLD, &req);
                }

                // Eliminate the for the local rows
                for (int elim_line = local_line + 1; elim_line < num_lines; elim_line++) {
                    // Get the scaling factor for elimination
                    double factor = local_matrix[elim_line][lines];

                    // Remove the pivot
                    for (int col = lines; col < N + 1; col++) {
                        local_matrix[elim_line][col] -=
                            local_matrix[local_line][col] * factor;
                    }
                }

                for (int i = check_rank + 1; i < procs; i++) {
                    MPI_Wait(&req, MPI_STATUS_IGNORE);
                }

            } else {
                // Receive pivot row
                MPI_Recv(pivot_line, N, MPI_DOUBLE, check_rank, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

                // Skip rows that have been fully processed
                for (int elim_row = 0; elim_row < num_lines; elim_row++) {
                    // Get the scaling factor for elimination
                    double factor = local_matrix[elim_row][lines];

                    // Remove the pivot
                    for (int col = lines; col < N + 1; col++) {
                        local_matrix[elim_row][col] -= pivot_line[col] * factor;
                    }
                }
            }
        }



        // recebe os resultados de cada sub matriz
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
