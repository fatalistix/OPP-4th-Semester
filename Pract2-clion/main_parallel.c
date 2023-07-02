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

void printIntVector(int * vector, int size)
{
    printf("[");
    for (int i = 0; i < size; ++i)
    {
        printf("%5d, ", vector[i]);
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

double countScalarSquare(double * vector, int size)
{
    double res = 0;
    for (int i = 0; i < size; ++i)
    {
        res += vector[i] * vector[i];
    }
    return res;
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
//    srand(time(NULL));
//    printf("time = %ld\n", time(NULL));
//    srand(1678250996);
    srand(1678536002);


    int matrixSize = 3000;
    int mpiRank;
    int mpiSize;
    int iterationsCounter = 0;


    double tao = 0.000001;
    double eps = 1e-5;
    double newEps;

    double iterationVectorScalarSquarePart = 0.;
    double iterationVectorScalarSquare = 0.;
    double vectorBScalarSquarePart = 0.;
    double vectorBScalarSquare = 0.;


    double * x;

    double * partOfA;
    double * partOfB;
    double * partOfItVec;

    int * parts;
    int * positions;

    const int rootRank = 0;




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
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if (mpiRank == rootRank)
    {
        double * A = (double *) malloc(sizeof(double) * matrixSize * matrixSize);
        double * b = (double *) malloc(sizeof(double) * matrixSize);

        x = (double *) malloc(sizeof(double) * matrixSize);

        makeMatrixRandomSymmetrical(A, matrixSize);
        makeVectorRandom(b, matrixSize, 100);
        makeVectorRandom(x, matrixSize, 10);

//        printMatrix(A, matrixSize, matrixSize);

        double start = MPI_Wtime();

        int minSegmentSize = matrixSize / mpiSize;
        int remains = matrixSize % mpiSize;

        // first will be filled for matrix
        parts     = (int *) malloc(sizeof(int) * mpiSize);
        positions = (int *) malloc(sizeof(int) * mpiSize);

        // filling parts and positions for matrix
        for (int i = 0; i < mpiSize; ++i)
        {
            parts[i]     = minSegmentSize * matrixSize;
            positions[i] = minSegmentSize * matrixSize * i;
            if (i < remains)
            {
                parts[i] += matrixSize;
                positions[i] += i * matrixSize;
            }
            else
            {
                positions[i] += remains * matrixSize;
            }
        }

        partOfA     = (double *) malloc(sizeof(double) * parts[mpiRank]);

        MPI_Scatterv(A, parts, positions, MPI_DOUBLE, partOfA, parts[mpiRank], MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        // filling parts and positions for vector
        for (int i = 0; i < mpiSize; ++i)
        {
            parts[i]     = minSegmentSize;
            positions[i] = minSegmentSize * i;
            if (i < remains)
            {
                parts[i]++;
                positions[i] += i;
            }
            else
            {
                positions[i] += remains;
            }
        }

//        printf("parts: ");
//        printIntVector(parts, mpiSize);

        partOfB     = (double *) malloc(sizeof(double) * parts[mpiRank]);
        partOfItVec = (double *) malloc(sizeof(double) * parts[mpiRank]);

        MPI_Scatterv(b, parts, positions, MPI_DOUBLE, partOfB, parts[mpiRank], MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
        vectorBScalarSquarePart = countScalarSquare(partOfB, parts[mpiRank]);
        MPI_Allreduce(&vectorBScalarSquarePart, &vectorBScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        newEps = eps * eps * vectorBScalarSquare;

//        printf("\n\n\nvb = %lf\n\n\n", vectorBScalarSquare);


        mulMatrixWithVector(partOfA, parts[mpiRank], matrixSize, x, partOfItVec);
        subVectors(partOfItVec, partOfB, partOfItVec, parts[mpiRank]);
//        printVector(partOfItVec, parts[mpiRank]);


        iterationVectorScalarSquarePart = countScalarSquare(partOfItVec, parts[mpiRank]);
        MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//        printf("===========================================\n");
//        printMatrix(partOfA, parts[mpiRank], matrixSize);
//        printf("===========================================\n");

//        printf("Positions: ");
//        printIntVector(positions, mpiSize);

//        printSLE(A, matrixSize, matrixSize, b);
//        printVector(x, matrixSize);



        while (iterationsCounter < 10000 && iterationVectorScalarSquare >= newEps)
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

        printf("Answer found for %d iterations\n Took %lf seconds\n", iterationsCounter, MPI_Wtime() - start);

        printf("b: ");
        printVector(b, matrixSize);

        mulMatrixWithVector(A, matrixSize, matrixSize, x, b);

        printf("Ax: ");
        printVector(b, matrixSize);


        free(A);
        free(b);
        free(x);
        free(parts);
        free(positions);
        free(partOfA);
        free(partOfB);
        free(partOfItVec);
    }



    else
    {
        x = (double *) malloc(sizeof(double) * matrixSize);

        int minSegmentSize = matrixSize / mpiSize;
        int remains = matrixSize % mpiSize;

        parts     = (int *) malloc(sizeof(int) * mpiSize);
        positions = (int *) malloc(sizeof(int) * mpiSize);

        for (int i = 0; i < mpiSize; ++i)
        {
            parts[i]     = minSegmentSize;
            positions[i] = minSegmentSize * i;
            if (i < remains)
            {
                parts[i]++;
                positions[i] += i;
            }
            else
            {
                positions[i] += remains;
            }
        }

        partOfA     = (double *) malloc(sizeof(double) * parts[mpiRank] * matrixSize);
        partOfB     = (double *) malloc(sizeof(double) * parts[mpiRank]);
        partOfItVec = (double *) malloc(sizeof(double) * parts[mpiRank]);

        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, partOfA, parts[mpiRank] * matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, partOfB, parts[mpiRank], MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
        vectorBScalarSquarePart = countScalarSquare(partOfB, parts[mpiRank]);
        MPI_Allreduce(&vectorBScalarSquarePart, &vectorBScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        newEps = eps * eps * vectorBScalarSquare;

        mulMatrixWithVector(partOfA, parts[mpiRank], matrixSize, x, partOfItVec);
        subVectors(partOfItVec, partOfB, partOfItVec, parts[mpiRank]);



        iterationVectorScalarSquarePart = countScalarSquare(partOfItVec, parts[mpiRank]);
        MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


        while (iterationsCounter < 10000 && iterationVectorScalarSquare >= newEps)
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

//        printf("===========================================\n");
//        printMatrix(partOfA, parts[mpiRank], matrixSize);
//        printf("===========================================\n");

        free(x);
        free(parts);
        free(positions);
        free(partOfA);
        free(partOfB);
        free(partOfItVec);
    }



    MPI_Finalize();
    return 0;
}
