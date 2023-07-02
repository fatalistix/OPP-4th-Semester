#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <omp.h>

#define SCHEDULE_MODE guided, 100000

void makeMatrixRandomSymmetrical(double * matrix, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = i; j < size; ++j)
        {
            matrix[i * size + j] = matrix[j * size + i] = 1. * ((rand() % 2) ? 1 : -1) * (rand() % 200) / (rand() % 200 + 1);               //1. * ((rand() % 2) ? 1 : -1) * rand() % 10 / (rand() % 10 + 1);
            if (i == j)
            {
                matrix[i * size + j] += 1.1 * size;
            }
        }
    }
}

void makeVectorZero(double * vector, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        vector[i] = 0.;
    }
}



void makeVectorRandom(double * vector, const int size, const int limit)
{
    for (int i = 0; i < size; ++i)
    {
        vector[i] = rand() % limit;
    }
}

void printSLE(const double * matrix, const int matrixHeight, const int matrixWidth, const double * vector)
{
    for (int i = 0; i < matrixHeight; ++i)
    {
        for (int j = 0; j < matrixWidth; ++j)
        {
            printf("%6.2lf ", matrix[i * matrixWidth + j]);
        }
        printf("| %6.2lf\n", vector[i]);
    }
}

void printMatrix(const double * matrix, const int matrixHeight, const int matrixWidth)
{
    for (int i = 0; i < matrixHeight; ++i)
    {
        for (int j = 0; j < matrixWidth; ++j)
        {
            printf("%6.2lf ", matrix[i * matrixWidth + j]);
        }
        printf("\n");
    }
}

void printVector(const double * vector, const int size)
{
    printf("[");
    for (int i = 0; i < size; ++i)
    {
        printf("%6.2lf, ", vector[i]);
    }
    printf("]\n");
}

void mulMatrixWithVector(const double * matrix, const int matrixHeight, const int matrixWidth, const double * vector, double * result)
{
    #pragma omp for schedule(SCHEDULE_MODE)
    for (int i = 0; i < matrixHeight; ++i)
    {
        result[i] = 0;
    }

    #pragma omp for /* collapse(2) */ schedule(SCHEDULE_MODE)
    for (int i = 0; i < matrixHeight; ++i)
    {
        for (int j = 0; j < matrixWidth; ++j)
        {
            // #pragma omp atomic
            result[i] += matrix[i * matrixWidth + j] * vector[j];
        }
    }
}

void sumVectors(const double * v1, const double * v2, double * res, const int size)
{
    #pragma omp for
    for (int i = 0; i < size; ++i)
    {
        res[i] = v1[i] + v2[i];
    }
}

void subVectors(const double * v1, const double * v2, double * res, const int size)
{
    #pragma omp for schedule(SCHEDULE_MODE)
    for (int i = 0; i < size; ++i)
    {
        res[i] = v1[i] - v2[i];
    }
}

void mulVectorWithScalar(const double * vector, double * res, const int size, const double scalar)
{
    #pragma omp for schedule(SCHEDULE_MODE)
    for (int i = 0; i < size; ++i)
    {
        res[i] = vector[i] * scalar;
    }
}



int simpleIterationMethod(const double * matrix, const double * vector, double * result, const int size, const double eps, const double tao)
{
    double * iterationVector = (double *) malloc(sizeof(double) * size); //* Ax^n - b

    double vectorBScalarSquare;
    double iterationVectorScalarSquare;
    double newEps;
    
    int iterationsCounter = 0;

    #pragma omp parallel
    {
        #pragma omp for reduction(+ : vectorBScalarSquare) schedule(SCHEDULE_MODE)
        for (int i = 0; i < size; ++i)
        {
            vectorBScalarSquare += vector[i] * vector[i];
        }

        mulMatrixWithVector(matrix, size, size, result, iterationVector);
        subVectors(iterationVector, vector, iterationVector, size);

        #pragma omp for reduction(+ : iterationVectorScalarSquare) schedule(SCHEDULE_MODE)
        for (int i = 0; i < size; ++i)
        {
            iterationVectorScalarSquare += iterationVector[i] * iterationVector[i];
        }

        #pragma omp single
        newEps = eps * eps * vectorBScalarSquare;

        #pragma omp barrier

        while (iterationsCounter < 10000 && iterationVectorScalarSquare >= newEps)
        {
            mulVectorWithScalar(iterationVector, iterationVector, size, tao);

            subVectors(result, iterationVector, result, size);

            mulMatrixWithVector(matrix, size, size, result, iterationVector);

            subVectors(iterationVector, vector, iterationVector, size);

            #pragma omp single
            ++iterationsCounter;

            #pragma omp single 
            iterationVectorScalarSquare = 0.;

            #pragma omp barrier

            #pragma omp for reduction(+ : iterationVectorScalarSquare) schedule(SCHEDULE_MODE)
            for (int i = 0; i < size; ++i)
            {
                iterationVectorScalarSquare += iterationVector[i] * iterationVector[i];
            }
        }
    }

    free(iterationVector);

    return iterationsCounter;
}


int main(int argc, char * argv[])
{
//    srand(time(NULL));
//    printf("time = %ld\n", time(NULL));
//    srand(1678250996);
    srand(1678536002);

    int matrixSize = 3000;
    double tao = 0.000001;
    int maxThreads = 1;
    if (argc > 1)
    {
        maxThreads = atoi(argv[1]);
    }
    if (argc > 2)
    {
        matrixSize = atoi(argv[2]);
    }
    if (argc > 3)
    {
        sscanf(argv[3], "%lf", &tao);
    }

    
    omp_set_num_threads(maxThreads);


    printf("\n\ntao = %lf\n\n", tao);



    double start = omp_get_wtime();

    double * A = (double *) malloc(sizeof(double) * matrixSize * matrixSize);
    double * b = (double *) malloc(sizeof(double) * matrixSize);
    double * x = (double *) malloc(sizeof(double) * matrixSize);

    makeMatrixRandomSymmetrical(A, matrixSize);
    makeVectorRandom(b, matrixSize, 100);
    makeVectorRandom(x, matrixSize, 10);

    // makeVectorZero(x, matrixSize);
    // printf("x: ");
    // printVector(x, matrixSize);

    printf("Answer found for %d iterations in %lf seconds\n", simpleIterationMethod(A, b, x, matrixSize, 1e-5, tao), omp_get_wtime() - start);
    // printf("Vector x: ");
    // printVector(x, matrixSize);

    // printf("SLE:\n");
    // printSLE(A, matrixSize, matrixSize, b);

    // double * buf = (double *) malloc(sizeof(double) * matrixSize);

    // mulMatrixWithVector(A, matrixSize, matrixSize, x, buf);
    // printf("Ax = ");
    // printVector(b, matrixSize);

    // int temp = findMismatchInVectors(b, buf, matrixSize);

    // if (temp == -1) { printf("\n\nVectors are equal\n\n"); }
    // else          { printf("\n\nVectors are different in %d pos\n\n", temp); }

    // printf("b: ");
    // printVector(b, matrixSize);
    // printf("Ax: ");
    // printVector(buf, matrixSize);


    // double * A = (double *) malloc(sizeof(double) * 9);
    // double * b = (double *) malloc(sizeof(double) * 3);
    // double * x = (double *) malloc(sizeof(double) * 3);

    // A[0] = 2; A[1] = 1; A[2] = 1;
    // A[3] = 1; A[4] = 2; A[5] = 1;
    // A[6] = 1; A[7] = 1; A[8] = 2;


    // b[0] = 4;
    // b[1] = 4;
    // b[2] = 4;


    free(A);
    free(b);
    free(x);
    //? free(buf);
    return 0;
}




// #!/bin/sh

// #PBS -l walltime=00:10:00
// #PBS -l select=1:ncpus=4:ompthreads=4
// #PBS -o out_parallel_new_4.txt
// #PBS -e err_parallel_new_4.txt

// cd $PBS_O_WORKDIR
// echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
// echo
// time ./parallel_new
