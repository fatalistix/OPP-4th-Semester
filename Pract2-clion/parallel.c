#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <mpi.h>

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

void printIntVector(int * vector, int size)
{
    printf("[");
    for (int i = 0; i < size; ++i)
    {
        printf("%5d, ", vector[i]);
    }
    printf("]\n");
}



int main(int argc, char * argv[])
{
    srand(1678536002);

    int matrixSize = 3000;
//    int matrixSize = 10;
    double tao = 0.000001;
    double eps = 1e-5;
    double newEps;
    double start;
    int minSegmentSize;
    int iterationsCounter = 0;

    double iterationVectorScalarSquarePart = 0.;
    double iterationVectorScalarSquare = 0.;
    double vectorBScalarSquarePart = 0.;
    double vectorBScalarSquare = 0.;


    const int rootRank = 0;

    double * A;
    double * b;
    double * x;

    double * partOfA;
    double * partOfB;
    double * partOfItVec;

    int * parts;
    int * positions;

    int mpiSize;
    int mpiRank;

    if (argc > 1)
    {
        matrixSize = atoi(argv[1]);
    }
    if (argc > 2)
    {
        sscanf(argv[2], "%lf", &tao);
    }

    printf("\n\ntao = %lf\n\n", tao);

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    x = (double *) malloc(sizeof(double) * matrixSize);
    parts = (int *) malloc(sizeof(int) * mpiSize);
    positions = (int *) malloc(sizeof(int) * mpiSize);

    if (mpiRank == rootRank)
    {
        A = (double *) malloc(sizeof(double) * matrixSize * matrixSize);
        b = (double *) malloc(sizeof(double) * matrixSize);

        makeMatrixRandomSymmetrical(A, matrixSize);
        makeVectorRandom(b, matrixSize, 100);
        makeVectorRandom(x, matrixSize, 10);



        parts[0] = matrixSize / mpiSize * matrixSize;
        positions[0] = 0;

        for (int i = 1; i < mpiSize; i++)
        {
            parts[i] = (matrixSize * (i + 1) / mpiSize - matrixSize * i / mpiSize) * matrixSize;
            positions[i] = parts[i - 1] + positions[i - 1];
        }
    }

    int segmentSize = matrixSize * (mpiRank + 1) / mpiSize - matrixSize * mpiRank / mpiSize;
    partOfA = (double *) malloc(sizeof(double) * segmentSize * matrixSize);

    MPI_Scatterv(A, parts, positions, MPI_DOUBLE, partOfA, segmentSize * matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);


    parts[0] = matrixSize / mpiSize;
    positions[0] = 0;

    for (int i = 1; i < mpiSize; i++)
    {
        parts[i] = (matrixSize * (i + 1) / mpiSize - matrixSize * i / mpiSize);
        positions[i] = parts[i - 1] + positions[i - 1];
    }

    if (mpiRank == rootRank)
    {
        printf("parts: ");
        printIntVector(parts, mpiSize);
        printf("positions: ");
        printIntVector(positions, mpiSize);
    }


    partOfB     = (double *) malloc(sizeof(double) * parts[mpiRank]);
    partOfItVec = (double *) malloc(sizeof(double) * parts[mpiRank]);

    MPI_Scatterv(b, parts, positions, MPI_DOUBLE, partOfB, parts[mpiRank], MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
    vectorBScalarSquarePart = countScalarSquare(partOfB, parts[mpiRank]);
    MPI_Allreduce(&vectorBScalarSquarePart, &vectorBScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    newEps = eps * eps * vectorBScalarSquare;

    mulMatrixWithVector(partOfA, parts[mpiRank], matrixSize, x, partOfItVec);
    subVectors(partOfItVec, partOfB, partOfItVec, parts[mpiRank]);

    iterationVectorScalarSquarePart = countScalarSquare(partOfItVec, parts[mpiRank]);
    MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    while (iterationsCounter < 100000 && iterationVectorScalarSquare >= newEps)
    {
        mulVectorWithScalar(partOfItVec, partOfItVec, parts[mpiRank], tao);
        subVectors(x + positions[mpiRank], partOfItVec, partOfItVec, parts[mpiRank]);

        MPI_Allgatherv(partOfItVec, parts[mpiRank], MPI_DOUBLE, x, parts, positions, MPI_DOUBLE, MPI_COMM_WORLD);

        mulMatrixWithVector(partOfA, parts[mpiRank], matrixSize, x, partOfItVec);
        subVectors(partOfItVec, partOfB, partOfItVec, parts[mpiRank]);

        iterationVectorScalarSquarePart = countScalarSquare(partOfItVec, parts[mpiRank]);
        MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        ++iterationsCounter;
    }

    if (mpiRank == rootRank)
    {
        printf("iter = %d, time = %lf\n", iterationsCounter, MPI_Wtime() - start);
        printf("b:  ");
        printVector(b, matrixSize);
        makeVectorRandom(b, matrixSize, 1);
        mulMatrixWithVector(A, matrixSize, matrixSize, x, b);
        printf("Ax: ");
        printVector(b, matrixSize);
    }

    MPI_Finalize();
    return 0;
}
// 6839