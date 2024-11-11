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

void gaussian_elimination(int n, double **m, int rank, int size) {
    double* pivo_linha;
    double** matriz;
    if(rank ==0){
        //bcast inicio
        for (size_t i = 0; i < n; i++)
        {
            double factor = m[i][i];
            m[i][i] = 1.0;

            for (int j = i + 1; j < n; j++) m[i][j] /= factor;
            MPI_Bcast(m[i], n , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        
    }
    else{

        MPI_Bcast(pivo_linha, n , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        int parcela = (n/size) * rank;
        int final = (rank == size-1)? n : parcela +(n/size); //ultimo processo pega oq sobrar
        for (int i = parcela; i < final; i++) { //linhas as quais aquele processo é responsável

            if (check_matrix(m, i)) i++;
            for (int j = i+1; j < n; j++)
            {
                double o = m[i][i];
                matriz[i][j] -= o* pivo_linha[j];
            }
            
        }

    }
for (int i = 0; i < n; i++) {

        if (check_matrix(m, i)) i++;


        double factor = m[i][i];
        m[i][i] = 1.0;

        for (int j = i + 1; j < n; j++) m[i][j] /= factor;

        //paralelizar
        for (int j = i + 1; j < n; j++) {
            double o = m[j][i];
            for (int k = 0; k < n; k++) m[j][k] = m[j][k] - o * m[i][k];
        }
    

}


int main(int argc, char *argv[]) {

    int meu_ranque, num_procs, inicio, dest, raiz=0;
    MPI_Status estado;

    srand((unsigned int) time(NULL));

    double **m = fill_matrix(N);

    order_matrix(N, m);

    gaussian_elimination(N, m);

    print_matrix(N, m);

    for (int i = 0; i < N; i++) free(m[i]);
    free(m);

    return 0;
}
