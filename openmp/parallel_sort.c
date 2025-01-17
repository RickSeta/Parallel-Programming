#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#define N 100000 // 5000 50000 100000

int *create_array(int const n) {
    int *m = malloc(sizeof(int) * n);

    for (int i = 0; i < n; i++) *(m + i) = rand() % (n + n) - n;

    return m;
}


void create_bitonic_seq(int *m, int n) {
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
    int *bubble_sort_seq = create_array(N);

    create_bitonic_seq(bitonic_seq, N);

    bubble_sort(bubble_sort_seq, N);

    free(bitonic_seq);
    free(bubble_sort_seq);

    return 0;
}
