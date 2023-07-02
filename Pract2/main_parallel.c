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
                matrix[i * size + j] += size * size;
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
    // srand(time(NULL));
    // srand(1677830020);
    srand(1678536002);


    int matrixSize = 10;
    if (argc > 1)
    {
        matrixSize = atoi(argv[1]);
    }

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

    //* This is matrix's height, but it will be also used for vectors
    int segmentSize;
    
    int iterationsCounter = 0;

    double iterationVectorScalarSquarePart = 0.;
    double iterationVectorScalarSquare = 0.;


    double vectorBScalarSquarePart = 0.;
    double vectorBScalarSquare = 0.;

    double eps = 1e-5;
    double tao = 0.01;
    double newEps;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    segmentSize  = (matrixSize / mpiSize); 

    //? Special rank that manages all memory and counts lesser then other
    int rootRank = mpiSize - 1;

    if (mpiRank == rootRank)
    {
        //? Usefull constant
        const int rootRankMatrixHeight = matrixSize - segmentSize * (mpiSize - 1);

        //*-------------------*//
        //* MEMORY ALLOCATION *//
        //*-------------------*//
        //* Allocating memory for matrix A, vector b and vector x
        //* Is used for init, for sending to other processes and for answer
        A = (double *) malloc(sizeof(double) * matrixSize * matrixSize);
        b = (double *) malloc(sizeof(double) * matrixSize);
        x = (double *) malloc(sizeof(double) * matrixSize);
        //* Manages separating parts for MPI_Scatterv func for matrix
        partsMatrix     = (int *) malloc(sizeof(int) * mpiSize);
        positionsMatrix = (int *) malloc(sizeof(int) * mpiSize);
        //* Manages separating parts for MPI_Scatterv func for vectors
        partsVector     = (int *) malloc(sizeof(int) * mpiSize);
        positionsVector = (int *) malloc(sizeof(int) * mpiSize);
        //* Allocating memory for part of matrix A, vector b and 
        //* iterationVector (all are used in calculations)
        partOfA               = (double *) malloc(sizeof(double) * rootRankMatrixHeight * matrixSize);
        partOfB               = (double *) malloc(sizeof(double) * rootRankMatrixHeight);
        partOfX               = (double *) malloc(sizeof(double) * rootRankMatrixHeight);
        partOfIterationVector = (double *) malloc(sizeof(double) * rootRankMatrixHeight);
        //*----------------------*//
        //* MEMORY ALLOCATION END*//
        //*----------------------*//

        //*----------------------------------*//
        //* ~~~~~~~~~~FILLING DATA~~~~~~~~~~ *//
        //*----------------------------------*//
        //* Sending separation's managing vectors
        for (int i = 0; i < mpiSize - 1; ++i)
        {
            partsMatrix    [i] = segmentSize * matrixSize;
            positionsMatrix[i] = segmentSize * matrixSize * i;
            partsVector    [i] = segmentSize;
            positionsVector[i] = segmentSize * i;
        }
        partsMatrix    [mpiSize - 1] = matrixSize * rootRankMatrixHeight;
        positionsMatrix[mpiSize - 1] = segmentSize * matrixSize * (mpiSize - 1);
        partsVector    [mpiSize - 1] = rootRankMatrixHeight;
        positionsVector[mpiSize - 1] = segmentSize * (mpiSize - 1);


        //* Initializing part of matrix A
        makeMatrixRandomSymmetrical(A, matrixSize);
        // fillWithZero(A + matrixSize * matrixSize, fictiousSize - matrixSize, matrixSize);
        
        //* Initializing vector b and x (answer vector)
        makeVectorRandom(b, matrixSize, 50);
        makeVectorRandom(x, matrixSize, 10);
        //*--------------------------------------*//
        //* ~~~~~~~~~~FILLING DATA END~~~~~~~~~~ *//
        //*--------------------------------------*//

        // printf("==============================\n");
        // printSLE(A, matrixSize, matrixSize, b);
        // printf("==============================\n");
        // printVector(x, matrixSize);
        // printf("==============================\n");
        // printf("rootRankMatrixHeight = %d; segmentSize = %d\n", rootRankMatrixHeight, segmentSize);
        // printf("==============================\n");
        // printf("partsMatrix: ");
        // printVectorInt(partsMatrix, mpiSize);
        // printf("positionsMatrix: ");
        // printVectorInt(positionsMatrix, mpiSize);
        // printf("partsVector: ");
        // printVectorInt(partsVector, mpiSize);
        // printf("positionsVector: ");
        // printVectorInt(positionsVector, mpiSize);
        // printf("==============================\n");

        //*----------------------------------*//
        //* ~~~~~~~~~~Scattering b~~~~~~~~~~ *//
        //*----------------------------------*//
        MPI_Scatterv(b, partsVector, positionsVector, MPI_DOUBLE, partOfB, rootRankMatrixHeight, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        // printf("==============================\n");
        // printf("b: ");
        // printVector(partOfB, rootRankMatrixHeight);
        // printf("==============================\n");

        //*----------------------------------------------*//
        //* ~~~~~~~~~~Counting b scalar square~~~~~~~~~~ *//
        //*----------------------------------------------*//
        vectorBScalarSquarePart = countScalarSquare(partOfB, rootRankMatrixHeight);
        // printf("==============================\n");
        // printf("Scalar square of part b = %lf\n", vectorBScalarSquarePart);
        // printf("==============================\n");
        MPI_Allreduce(&vectorBScalarSquarePart, &vectorBScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // printf("==============================\n");
        // printf("Scalar square of part b = %lf, Scalar square of b = %lf\n", vectorBScalarSquarePart, vectorBScalarSquare);
        // printf("==============================\n");

        //*-------------------------------------*//
        //* ~~~~~~~~~~Counting newEps~~~~~~~~~~ *//
        //*-------------------------------------*//
        newEps = eps * eps * vectorBScalarSquare;

        //*----------------------------------*//
        //* ~~~~~~~~~~Scattering A~~~~~~~~~~ *//
        //*----------------------------------*//
        MPI_Scatterv(A, partsMatrix, positionsMatrix, MPI_DOUBLE, partOfA, rootRankMatrixHeight * matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        printf("==============================\n");
        printMatrix(partOfA, rootRankMatrixHeight, matrixSize);
        printf("==============================\n");

        //*------------------------------------*//
        //* ~~~~~~~~~~Broadcasting x~~~~~~~~~~ *//
        //*------------------------------------*//
        MPI_Bcast(x, matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        //*----------------------------------------------*//
        //* ~~~~~~~~~~Counting IterationVector~~~~~~~~~~ *//
        //*----------------------------------------------*//
        // printf("==============================\n");
        // printf("x: ");
        // printVector(x, matrixSize);
        // printf("==============================\n");
        mulMatrixWithVector(partOfA, rootRankMatrixHeight, matrixSize, x, partOfIterationVector);
        // printf("==============================\n");
        // printf("(partA)x: ");
        // printVector(partOfIterationVector, rootRankMatrixHeight);
        // printf("==============================\n");
        subVectors(partOfIterationVector, partOfB, partOfIterationVector, rootRankMatrixHeight);
        // printf("==============================\n");
        // printf("(partA)x-(partB): ");
        // printVector(partOfIterationVector, rootRankMatrixHeight);
        // printf("==============================\n");
        //*-------------------------------------------------------------*//
        //* ~~~~~~~~~~Counting iterationVector's scalarSquare~~~~~~~~~~ *//
        //*-------------------------------------------------------------*//
        iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, rootRankMatrixHeight);
        // printf("==============================\n");
        // printf("itVecScalPart = %lf\n", iterationVectorScalarSquarePart);
        // printf("==============================\n");
        MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        while (iterationsCounter < 50 && iterationVectorScalarSquare >= newEps)
        {
            // printf("%lf on %d on iteration %d\n", iterationVectorScalarSquare, mpiRank, iterationsCounter);
            mulVectorWithScalar(partOfIterationVector, partOfIterationVector, rootRankMatrixHeight, tao);
            subVectors(x + positionsVector[mpiRank], partOfIterationVector, partOfX, rootRankMatrixHeight);

            MPI_Allgatherv(partOfX, rootRankMatrixHeight, MPI_DOUBLE, x, partsVector, positionsVector, MPI_DOUBLE, MPI_COMM_WORLD);

            mulMatrixWithVector(partOfA, rootRankMatrixHeight, matrixSize, x, partOfIterationVector);
            subVectors(partOfIterationVector, partOfB, partOfIterationVector, rootRankMatrixHeight);

            iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, rootRankMatrixHeight);
            MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            ++iterationsCounter;
        }


        printf("Answer found for %d iterations\n", iterationsCounter);




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
        x = (double *) malloc(sizeof(double) * matrixSize);
        partOfA = (double *) malloc(sizeof(double) * matrixSize * segmentSize);
        partOfB = (double *) malloc(sizeof(double) * segmentSize);
        partOfX = (double *) malloc(sizeof(double) * segmentSize);
        partOfIterationVector = (double *) malloc(sizeof(double) * segmentSize);

        //? better to do it without scattering, anyway it will ne initialized on zero so we will parallelly do it
        int * partsVector     = (int *) malloc(sizeof(int) * mpiSize);
        int * positionsVector = (int *) malloc(sizeof(int) * mpiSize);

        for (int i = 0; i < mpiSize - 1; ++i)
        {
            partsVector    [i] = segmentSize;
            positionsVector[i] = segmentSize * i;
        }
        partsVector    [mpiSize - 1] = matrixSize - segmentSize * (mpiSize - 1);
        positionsVector[mpiSize - 1] = segmentSize * (mpiSize - 1);



        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, partOfB, segmentSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        vectorBScalarSquarePart = countScalarSquare(partOfB, segmentSize);
        MPI_Allreduce(&vectorBScalarSquarePart, &vectorBScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        newEps = eps * eps * vectorBScalarSquare;

        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, partOfA, segmentSize * matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        MPI_Bcast(x, matrixSize, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

        mulMatrixWithVector(partOfA, segmentSize, matrixSize, x, partOfIterationVector);
        subVectors(partOfIterationVector, partOfB, partOfIterationVector, segmentSize);

        iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, segmentSize);
        MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // printf("iteration vector on %d: ", mpiRank);
        // printVector(partOfIterationVector, segmentSize);

        while (iterationsCounter < 50 && iterationVectorScalarSquare >= newEps)
        {
            // printf("%lf on %d on iteration %d\n", iterationVectorScalarSquare, mpiRank, iterationsCounter);
            mulVectorWithScalar(partOfIterationVector, partOfIterationVector, segmentSize, tao);
            subVectors(x + mpiRank * segmentSize, partOfIterationVector, partOfX, segmentSize);

            MPI_Allgatherv(partOfX, segmentSize, MPI_DOUBLE, x, partsVector, positionsVector, MPI_DOUBLE, MPI_COMM_WORLD);

            mulMatrixWithVector(partOfA, segmentSize, matrixSize, x, partOfIterationVector);
            subVectors(partOfIterationVector, partOfB, partOfIterationVector, segmentSize);

            iterationVectorScalarSquarePart = countScalarSquare(partOfIterationVector, segmentSize);
            MPI_Allreduce(&iterationVectorScalarSquarePart, &iterationVectorScalarSquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            ++iterationsCounter;
        }

    
        free(x);
        free(partOfA);
        free(partOfB);
        free(partOfX);
        // printf("]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\n");
        // printf("rank = %d\n", mpiRank);
        // printf("]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\n");
        free(partOfIterationVector);
        free(partsVector);
        free(positionsVector);
    }

    MPI_Finalize();
    return 0;
}
