#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

void makeMatrixRandomSymmetrical(double * matrix, int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = i; j < size; ++j)
        {
            matrix[i * size + j] = matrix[j * size + i] = 1. * ((rand() % 2) ? 1 : -1) * (rand() % 200) / (rand() % 200 + 1);               //1. * ((rand() % 2) ? 1 : -1) * rand() % 10 / (rand() % 10 + 1);
            if (i == j)
            {
                matrix[i * size + j] += size / 2;
            }
        }
    }
}

//void makeMatrixRandomSymmetrical(double * matrix, int size)
//{
//    for (int i = 0; i < size; ++i)
//    {
//        for (int j = i; j < size; ++j)
//        {
//            matrix[i * size + j] = matrix[j * size + i] = j % 1000;         // * ((rand() % 2) ? 1 : -1);
//            if (i == j)
//            {
//                matrix[i * size + j] += size;
//            }
//        }
//    }
//}

void makeVectorZero(double * vector, int size)
{
    for (int i = 0; i < size; ++i)
    {
        vector[i] = 0.;
    }
}



void makeVectorRandom(double * vector, int size, int limit)
{
    for (int i = 0; i < size; ++i)
    {
        vector[i] = rand() % limit;
    }
}

void printSLE(double * matrix, int matrixHeight, int matrixWidth, double * vector)
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

void printMatrix(double * matrix, int matrixHeight, int matrixWidth)
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

void printVector(double * vector, int size)
{
    printf("[");
    for (int i = 0; i < size; ++i)
    {
        printf("%6.2lf, ", vector[i]);
    }
    printf("]\n");
}

void mulMatrixWithVector(double * matrix, int matrixHeight, int matrixWidth, double * vector, double * result)
{
    for (int i = 0; i < matrixHeight; ++i)
    {
        result[i] = 0;
        for (int j = 0; j < matrixWidth; ++j)
        {
            result[i] += matrix[i * matrixWidth + j] * vector[j];
        }
    }
}

void sumVectors(double * v1, double * v2, double * res, int size)
{
    for (int i = 0; i < size; ++i)
    {
        res[i] = v1[i] + v2[i];
    }
}

void subVectors(double * v1, double * v2, double * res, int size)
{
    for (int i = 0; i < size; ++i)
    {
        res[i] = v1[i] - v2[i];
    }
}

void mulVectorWithScalar(double * vector, double * res, int size, double scalar)
{
    for (int i = 0; i < size; ++i)
    {
        res[i] = vector[i] * scalar;
    }
}

double countScalarSquare(double * vector, int size)
{
    double res = 0;
    for (int i = 0; i < size; ++i)
    {
        res += vector[i] * vector[i];
    }
    return res;
}



int simpleIterationMethod(double * matrix, double * vector, double * result, int size, double eps, double tao)
{
    printf("size=%d\n\n\n", size);

    int iterationsCounter = 0;
    double vectorBScalarSquare = countScalarSquare(vector, size);
    double * iterationVector = (double *) malloc(sizeof(double) * size); //* Ax^n - b

    mulMatrixWithVector(matrix, size, size, result, iterationVector);
    // printf("Ax: ");
    // printVector(iterationVector, size);
    subVectors(iterationVector, vector, iterationVector, size);

    // printf("Ax-b: ");
    // printVector(iterationVector, size);
    double newEps = eps * eps * vectorBScalarSquare;

    while (iterationsCounter < 5000 && countScalarSquare(iterationVector, size) >= newEps)
    {
        // printf("%lf on iteration %d\n", countScalarSquare(iterationVector, size), iterationsCounter);
        // printf("iteration %d:\n", iterationsCounter);
        // printf("iterationVector: ");
        // printVector(iterationVector, size);
        mulVectorWithScalar(iterationVector, iterationVector, size, tao);
        // printf("iterationVector after mulVectorWithScalar: ");
        // printVector(iterationVector, size);
        // printf("result vector: ");
        // printVector(result, size);
        subVectors(result, iterationVector, result, size);


        mulMatrixWithVector(matrix, size, size, result, iterationVector);
        // printf("iteration vector %d: ", iterationsCounter);
        // printVector(iterationVector, size);
        subVectors(iterationVector, vector, iterationVector, size);

        // printf("\n\n\nresult vector at %d iter end:\n", iterationsCounter);
        // printVector(result, size);
        // printf("\n\n\n");

        ++iterationsCounter;
    }
    free(iterationVector);
    return iterationsCounter;
}



int findMismatchInVectors(double * v1, double * v2, int size)
{
    for (int i = 0; i < size; ++i)
    {
        if (v1[i] - v2[i] >= 0.01)
        {
            return i;
        }
    }
    return -1;
}



int main(int argc, char * argv[])
{
    // srand(time(NULL));
    // printf("time = %ld\n", time(NULL));
    srand(1678250996);


    int matrixSize = 2500;
    if (argc > 1)
    {
        matrixSize = atoi(argv[1]);
    }



    double tao = 0.00001;
    if (argc > 2)
    {
        sscanf(argv[2], "%lf", &tao);
    }

    printf("\n\ntao = %lf\n\n", tao);



    double start;
    MPI_Init(&argc, &argv);
    start = MPI_Wtime();

    double * A = (double *) malloc(sizeof(double) * matrixSize * matrixSize);
    double * b = (double *) malloc(sizeof(double) * matrixSize);
    double * x = (double *) malloc(sizeof(double) * matrixSize);

    makeMatrixRandomSymmetrical(A, matrixSize);
    makeVectorRandom(b, matrixSize, 500);
    makeVectorRandom(x, matrixSize, 10);
    // makeVectorZero(x, matrixSize);
    // printf("x: ");
    // printVector(x, matrixSize);

    printf("Answer found for %d in iterations in %lf seconds\n", simpleIterationMethod(A, b, x, matrixSize, 1e-5, tao), MPI_Wtime() - start);
    // printf("Vector x: ");
    // printVector(x, matrixSize);

    // printf("SLE:\n");
    // printSLE(A, matrixSize, matrixSize, b);

    double * buf = (double *) malloc(sizeof(double) * matrixSize);

    mulMatrixWithVector(A, matrixSize, matrixSize, x, buf);
    // printf("Ax = ");
    // printVector(b, matrixSize);

    int temp = findMismatchInVectors(b, buf, matrixSize);

    if (temp == -1) { printf("\n\nVectors are equal\n\n"); }
    else          { printf("\n\nVectors are different in %d pos\n\n", temp); }

    printf("b: ");
    printVector(b, matrixSize);
    printf("Ax: ");
    printVector(buf, matrixSize);


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
    free(buf);

    MPI_Finalize();
    return 0;
}

