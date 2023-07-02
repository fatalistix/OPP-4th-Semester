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
                matrix[i * size + j] += 1.1 * size;
            }
        }
    }
}

void makeVectorRandom(double * vector, int size, int limit)
{
    for (int i = 0; i < size; ++i)
    {
        vector[i] = rand() % limit;
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

void printVectorInt(int * vector, int size)
{
    printf("[");
    for (int i = 0; i < size; ++i)
    {
        printf("%6d, ", vector[i]);
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

int main(int argc, char * argv[])
{
    srand(1678536002);

    int matrixSize = 3000;

    double * A = NULL;
    double * b = NULL;
    double * x = NULL;

    double * partOfA = NULL;
    double * partOfB = NULL;
    double * partOfX = NULL;
    double * partOfIterationVector = NULL;

    int * partsMatrix;
    int * positionsMatrix;
    int * partsVector;
    int * positionsVector;

    int mpiRank;
    int mpiSize;

    int iterationsCounter = 0;

    double iterationVectorScalarSquarePart = 0.;
    double iterationVectorScalarSquare = 0.;


    double vectorBScalarSquarePart = 0.;
    double vectorBScalarSquare = 0.;

    double eps = 1e-5;
    double tao = 0.000001;
    double newEps;


    int rootRank = 0;
    double start;


    if (argc > 1)
    {
        matrixSize = atoi(argv[1]);
    }
    if (argc > 2)
    {
        sscanf(argv[2], "%lf", &tao);
    }

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if (mpiRank == rootRank)
    {
        A = (double *) malloc(sizeof(double) * matrixSize * matrixSize);
        b = (double *) malloc(sizeof(double) * matrixSize);
        x = (double *) malloc(sizeof(double) * matrixSize);

        partsMatrix     = (int *) malloc(sizeof(int) * mpiSize);
        positionsMatrix = (int *) malloc(sizeof(int) * mpiSize);

        partsVector     = (int *) malloc(sizeof(int) * mpiSize);
        positionsVector = (int *) malloc(sizeof(int) * mpiSize);


        partsMatrix[0] = matrixSize / mpiSize * matrixSize;
        partsVector[0] = matrixSize / mpiSize;
        positionsMatrix[0] = 0;
        positionsVector[0] = 0;

        for (int i = 1; i < mpiSize; i++)
        {
            partsMatrix[i] = (matrixSize * (i + 1) / mpiSize - matrixSize * i / mpiSize) * matrixSize;
            partsVector[i] = (matrixSize * (i + 1) / mpiSize - matrixSize * i / mpiSize);
            positionsMatrix[i] = partsMatrix[i - 1] + positionsMatrix[i - 1];
            positionsVector[i] = partsVector[i - 1] + positionsVector[i - 1];
        }

        partOfA               = (double *) malloc(sizeof(double) * partsMatrix[mpiRank]);
        partOfB               = (double *) malloc(sizeof(double) * partsVector[mpiRank]);
        partOfX               = (double *) malloc(sizeof(double) * partsVector[mpiRank]);
        partOfIterationVector = (double *) malloc(sizeof(double) * partsVector[mpiRank]);

        makeMatrixRandomSymmetrical(A, matrixSize);
        makeVectorRandom(b, matrixSize, 100);
        makeVectorRandom(x, matrixSize, 10);

        MPI_Scatterv(b, partsVector, positionsVector, MPI_DOUBLE, partOfB, partsVector[mpiRank], MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        vectorBScalarSquarePart = countScalarSquare(partOfB, partsVector[mpiRank]);
        MPI_Allreduce(&vectorBScalarSquarePart, &vectorBScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        newEps = eps * eps * vectorBScalarSquare;

        MPI_Scatterv(A, partsMatrix, positionsMatrix, MPI_DOUBLE, partOfA, partsMatrix[mpiRank], MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        MPI_Bcast(x, matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        mulMatrixWithVector(partOfA, partsVector[mpiRank], matrixSize, x, partOfIterationVector);
        subVectors(partOfIterationVector, partOfB, partOfIterationVector, partsVector[mpiRank]);

        iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, partsVector[mpiRank]);
        MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        while (iterationsCounter < 10000 && iterationVectorScalarSquare >= newEps)
        {
            mulVectorWithScalar(partOfIterationVector, partOfIterationVector, partsVector[mpiRank], tao);
            subVectors(x + positionsVector[mpiRank], partOfIterationVector, partOfX, partsVector[mpiRank]);

            MPI_Allgatherv(partOfX, partsVector[mpiRank], MPI_DOUBLE, x, partsVector, positionsVector, MPI_DOUBLE, MPI_COMM_WORLD);

            mulMatrixWithVector(partOfA, partsVector[mpiRank], matrixSize, x, partOfIterationVector);
            subVectors(partOfIterationVector, partOfB, partOfIterationVector, partsVector[mpiRank]);

            iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, partsVector[mpiRank]);
            MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            ++iterationsCounter;
        }

        printf("Answer found for %d iterations\n", iterationsCounter);

        printf("b:  ");
        printVector(b, matrixSize);
        makeVectorRandom(b, matrixSize, 1);
        mulMatrixWithVector(A, matrixSize, matrixSize, x, b);
        printf("Ax: ");
        printVector(b, matrixSize);


        free(A);
        free(b);
        free(x);
        free(partsMatrix);
        free(positionsMatrix);
        free(partsVector);
        free(positionsVector);
        free(partOfA);
        free(partOfB);
        free(partOfX);
        free(partOfIterationVector);
    }

    else
    {
        partsVector     = (int *) malloc(sizeof(int) * mpiSize);
        positionsVector = (int *) malloc(sizeof(int) * mpiSize);

        partsVector[0] = matrixSize / mpiSize;
        positionsVector[0] = 0;

        for (int i = 1; i < mpiSize; i++)
        {
            partsVector[i] = (matrixSize * (i + 1) / mpiSize - matrixSize * i / mpiSize);
            positionsVector[i] = partsVector[i - 1] + positionsVector[i - 1];
        }

        x = (double *) malloc(sizeof(double) * matrixSize);
        partOfA = (double *) malloc(sizeof(double) * partsVector[mpiRank] * matrixSize);
        partOfB = (double *) malloc(sizeof(double) * partsVector[mpiRank]);
        partOfX = (double *) malloc(sizeof(double) * partsVector[mpiRank]);
        partOfIterationVector = (double *) malloc(sizeof(double) * partsVector[mpiRank]);

        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, partOfB, partsVector[mpiRank], MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        vectorBScalarSquarePart = countScalarSquare(partOfB, partsVector[mpiRank]);
        MPI_Allreduce(&vectorBScalarSquarePart, &vectorBScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        newEps = eps * eps * vectorBScalarSquare;

        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, partOfA, partsVector[mpiRank] * matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        MPI_Bcast(x, matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        mulMatrixWithVector(partOfA, partsVector[mpiRank], matrixSize, x, partOfIterationVector);
        subVectors(partOfIterationVector, partOfB, partOfIterationVector, partsVector[mpiRank]);

        iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, partsVector[mpiRank]);
        MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        while (iterationsCounter < 10000 && iterationVectorScalarSquare >= newEps)
        {
            mulVectorWithScalar(partOfIterationVector, partOfIterationVector, partsVector[mpiRank], tao);
            subVectors(x + positionsVector[mpiRank], partOfIterationVector, partOfX, partsVector[mpiRank]);

            MPI_Allgatherv(partOfX, partsVector[mpiRank], MPI_DOUBLE, x, partsVector, positionsVector, MPI_DOUBLE, MPI_COMM_WORLD);

            mulMatrixWithVector(partOfA, partsVector[mpiRank], matrixSize, x, partOfIterationVector);
            subVectors(partOfIterationVector, partOfB, partOfIterationVector, partsVector[mpiRank]);

            iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, partsVector[mpiRank]);
            MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            ++iterationsCounter;
        }

        free(x);
        free(partOfA);
        free(partOfB);
        free(partOfX);
        free(partOfIterationVector);
        free(partsVector);
        free(positionsVector);
    }

    MPI_Finalize();
    return 0;
}
