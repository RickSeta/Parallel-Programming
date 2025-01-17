#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
typedef enum { false, true } bool;

#define N 4 // 5000 50000 100000

int *create_array(int const n) {
    int *m = malloc(sizeof(int) * n);

    for (int i = 0; i < n; i++) *(m + i) = rand() % 1000;

    return m;
}

int* merge(int *arr, int menor_index, int tamanho, bool ordem_crescente) {
    if(tamanho > 1){
        int meio = tamanho / 2;

        for (int i = menor_index; i < menor_index + meio; i++) {
            if (ordem_crescente == (arr[i] > arr[i + meio])) {
                int temp = arr[i];
                arr[i] = arr[i + meio];
                arr[i + meio] = temp;
            }
        }
        
        merge(arr, menor_index, meio, ordem_crescente);
        merge(arr, menor_index + meio, meio, ordem_crescente);
    }
}

int* create_bitonic_seq(int *arr, int menor_index, int tamanho, bool ordem_crescente) {
    if(tamanho > 1){
        int meio = tamanho / 2;
        create_bitonic_seq(arr, menor_index, meio, true);
        create_bitonic_seq(arr, menor_index + meio, meio, false);

        merge(arr, menor_index, tamanho, ordem_crescente);
    }
}



void bubble_sort(int *m, const int n) {
    clock_t start = clock();

    for (int i = 0; i < n; i++) {
#pragma omp parallel for num_threads(2), default(none) shared(i,n,m)
        for (int j = i % 2; j < n - 1; j += 2) {
            if (m[j] > m[j + 1]) {
                int const temp = m[j];
                m[j] = m[j + 1];
                m[j + 1] = temp;
            }
        }
    }

    clock_t end = clock();

    printf("Tempo Gasto: %.4f s.\n", (double) (end - start) / (double) CLOCKS_PER_SEC);
}


int main() {
    srand(time(NULL));

    int *bitonic_seq = create_array(N);
    // int *bubble_sort_seq = create_array(N);
    for (int i = 0; i < N; i++) {
        printf("%d ", bitonic_seq[i]);
    }
    printf("\n");
    create_bitonic_seq(bitonic_seq, 0, N, true);

    // bubble_sort(bubble_sort_seq, N);
    for (int i = 0; i < N; i++) {
        printf("%d ", bitonic_seq[i]);
    }
    printf("\n");
    free(bitonic_seq);
    // free(bubble_sort_seq);

    return 0;
}