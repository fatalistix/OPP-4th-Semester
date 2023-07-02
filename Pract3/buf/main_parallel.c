#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define PLANNING_TYPE auto

static const double TAO = 0.0001;
static const double EPSILON = 0.0000001;

void initA(double *A, const int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            if (i != j) {
                A[i * N + j] = A[j * N + i] = ((double) rand() / RAND_MAX) * 2 - 1;
            } else {
                A[i * N + j] = N / 8.;
            }
        }
    }
}

void initB(double *b, const int N) {
    for (int i = 0; i < N; i++) {
        b[i] = N + 1;
    }
}

void initX(double *x, const int N) {
    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }
}

void
mulVectorOnScalar(const double *vector, double *scalarMatrix, const int N, const double scalar) {
#pragma omp for schedule(PLANNING_TYPE)
    for (int i = 0; i < N; i++) {
        scalarMatrix[i] = scalar * vector[i];
    }
}

void
subVectors(const double *vector1, const double *vector2, double *subVector, const int N) {
#pragma omp for schedule(PLANNING_TYPE)
    for (int i = 0; i < N; i++) {
        subVector[i] = vector1[i] - vector2[i];
    }
}

void mulMatrixOnVector(const double *matrix, const double *vector, double *mulVector, const int N) {
#pragma omp for schedule(PLANNING_TYPE)
    for (int i = 0; i < N; i++) {
        mulVector[i] = 0;
    }

#pragma omp for collapse(2) schedule(PLANNING_TYPE)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
#pragma omp atomic
            mulVector[i] += matrix[i * N + j] * vector[j];
        }
    }
}

double vectorModule(const double *vector, const int N) {
    double sum;

#pragma omp single
    {
        sum = 0;
    }

#pragma omp barrier

#pragma omp shared(sum) for reduction(+ : critVectModule) schedule(PLANNING_TYPE)
    for (int i = 0; i < N; i++) {
        sum += vector[i] * vector[i];
        printf("i am %d in countScalarSquare loop\n", omp_get_thread_num());
    }

#pragma omp barrier
    return sum;
}

int main(int argc, char **argv) {
    srand(1792);
    const int N = atoi(argv[1]);
    const int numThreads = atoi(argv[2]);

    double *A = (double *) malloc(N * N * sizeof(double));
    double *b = (double *) malloc(N * sizeof(double));
    double *xn = (double *) malloc(N * sizeof(double));
    double *x = (double *) malloc(N * sizeof(double));

    double end, start;

    initA(A, N);
    initB(b, N);
    initX(xn, N);

    start = omp_get_wtime();

    double critVectModule;
    double newEpsilon = vectorModule(b, N) * EPSILON * EPSILON;

    omp_set_num_threads(numThreads);

#pragma omp parallel
    {

        int countIterations = 0;

        do {

            mulMatrixOnVector(A, xn, x, N);

            subVectors(x, b, x, N);

#pragma omp barrier

            critVectModule = vectorModule(x, N);

            mulVectorOnScalar(x, x, N, TAO);

            subVectors(xn, x, x, N);

#pragma omp barrier

#pragma omp master
            {
                memcpy(xn, x, N * sizeof(double));
            }

#pragma omp barrier

            ++countIterations;

        } while (newEpsilon <= critVectModule && countIterations < 50000);

#pragma omp master
        {
            end = omp_get_wtime();
            printf("OMP_NUM_THREADS = 16 \n");
            printf("Iterations: %d\n", countIterations);
            printf("Time taken: %f sec \n", end - start);
        }
    }

    free(A);
    free(b);
    free(xn);
    free(x);

    return 0;
}
