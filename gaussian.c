#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
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

double** appendMatrix(int n, double matrix[n][n], double vector[n]) {
    int rows = n;
    int cols = n;
    int new_cols = cols + 1;

    // Allocate memory for the new matrix
    double** new_matrix = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        new_matrix[i] = (double*)malloc(new_cols * sizeof(double));
    }

    // Copy elements from the original matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            new_matrix[i][j] = matrix[i][j];
        }
    }

    // Copy elements from the vector to the new column
    for (int i = 0; i < rows; i++) {
        new_matrix[i][new_cols - 1] = vector[i];
    }

    return new_matrix;
}

void gaussianElimination(int n, double** matrix) {
    // Augment the matrix with the vector

    


    
}

int main(int argc, char **argv[]){

    int meu_ranque, num_procs, inicio, dest, raiz=0;

    int n = 3; // Number of variables
    double matrix[3][3] = {
        {2, 1, -1},
        {-3, -1, 2},
        {-2, 1, 2}
    };
    double vector[3] = {8, -11, -3};
    double** expanded_matrix = appendMatrix(3, matrix, vector);




    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n+1; j++) {
            printf("%lf ", expanded_matrix[i][j]);
        }
        printf("\n");
    }

    // Free the memory
    for (int i = 0; i < n; i++) {
        free(expanded_matrix[i]);
    }
    free(expanded_matrix);

    return 0;
}
