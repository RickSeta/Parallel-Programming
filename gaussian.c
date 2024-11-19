#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>
#define N 5

void print_matrix(int n, double **m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%8.2lf ", m[i][j]);
        }
        puts("");
    }
    puts("");
}

int check_matrix(double **m, int n) {
    return m[n][n] == 0.0;
}

double **fill_matrix(int n) {

    double **m = (double **) malloc(sizeof(double *) * n);

    for (int i = 0; i < n; i++)
        m[i] = (double *) malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            //m[i][j] = 3*(i==j);
            m[i][j] = rand() / (double) RAND_MAX * (10 + 8) - 8.0;
        }
    }

    return m;

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

    //paralelizar

    for (int i = 0; i < n-1; i++) {

        if (check_matrix(m, i)) i++;

        double max_value = m[i][i];
        int max_value_pos = i;
        for (j = n - 1; j > k; j--) {
            if (fabs(m[j][i]) > max_value) {
                max_value = m[j][i];
                max_value_pos = j;
            }
        }
        if (max_value != m[i][i])
            change_rows(n, m, max_value_pos, j);
        k++;
    }

}

void gaussian_elimination(int n_linhas, double **m, int rank, int num_proc) {
    double* pivo_linha;
    for (int firstLine = 0; firstLine < n_linhas; firstLine++) {

        
        if (rank == 0) {
            // Prepare pivot row in process 0
            pivo_linha = m[firstLine];
            double piv = pivo_linha[firstLine];

            for (int i = 0; i < n_linhas; i++){ 
                printf("%lf ",pivo_linha[i]);
                pivo_linha[i] /= piv;
                
                printf("%lf ",pivo_linha[i]);
                }
            
            printf("\n ");
        }
        MPI_Bcast(pivo_linha, n_linhas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //todo mundo tem a linha pivo a partir daqui

        // Dynamically calculate the workload for each process in this iteration
        int rows_left = n_linhas - firstLine - 1; // Remaining rows after the pivot row
        int rows_per_proc = rows_left / num_proc; // Basic share of rows per process
        int extra_rows = rows_left % num_proc; // Extra rows to distribute

        int start_row = firstLine + 1 + rank * rows_per_proc;
        int end_row = start_row + rows_per_proc + ((rank == num_proc-1)? extra_rows: 0);

        //cada processo tem o numero certo de linhas a partir daqui

        // Process each assigned row within this iteration of firstLine
        for (int i = start_row; i < end_row; i++) {
            //if (check_matrix(m, i)) continue;
            double factor = m[i][firstLine];
            for (int j = firstLine + 1; j < n_linhas; j++) {
                m[i][j] -= factor * pivo_linha[j];
            }
        }
    }

}

int main(int argc, char **argv[]){

    int meu_ranque, num_procs, inicio, dest, raiz=0;
    MPI_Status estado;

    srand((unsigned int) time(NULL));

    double **m = fill_matrix(N);

    order_matrix(N, m);

    MPI_Init(&argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &meu_ranque);

    gaussian_elimination(N, m, meu_ranque, num_procs);
    
    MPI_Finalize();


    for (int i = 0; i < N; i++) free(m[i]);
    free(m);

    return 0;
}
