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
            matrix[i * size + j] = matrix[j * size + i] = 1. * ((rand() % 2) ? 1 : -1) * (rand() % 200) / (rand() % 200 + 1);
            if (i == j)
            {
                matrix[i * size + j] += size / 2;
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
    // srand(time(NULL));
    // printf("time = %ld\n", time(NULL));
    srand(1678250996);


    int matrixSize = 2500;
    int mpiRank;
    int mpiSize;
    int iterationsCounter = 0;


    double tao = 0.00001;
    double eps = 1e-5;
    double newEps;

    double iterationVectorScalarSquarePart = 0.;
    double iterationVectorScalarSquare = 0.;
    double vectorBScalarSquarePart = 0.;
    double vectorBScalarSquare = 0.;


    double * x;

    double * partOfA;
    double * partOfB;
    double * partOfX; 
    double * partOfIterationVector;

    int * partsVector;
    int * positionsVector;




    if (argc > 1)
    {
        matrixSize = atoi(argc[1]);
    }
    if (argc > 2)
    {
        sscanf(argv[2], "%lf", &tao);
    }

    printf("\n\ntao = %lf\n\n", tao);
    
    
    MPI_Init(&argc, &argv[]);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if (mpiRank == 0)
    {
        double * A = (double *) malloc(sizeof(double) * matrixSize * matrixSize);
        double * b = (double *) malloc(sizeof(double) * matrixSize);

        x = (double *) malloc(sizeof(double) * matrixSize);

        makeMatrixRandomSymmetrical(A, matrixSize);
        makeVectorRandom(b, matrixSize, 500);
        makeVectorRandom(x, matrixSize, 3)

        double start = MPI_Wtime();

        int minSegmentSize = matrixSize / mpiSize;
        int remains = matrixSize % mpiSize;

        int * partsMatrix     = (int *) malloc(sizeof(int) * mpiSize);
        int * positionsMatrix = (int *) malloc(sizeof(int) * mpiSize);

        
    }


    else
    {

    }



    MPI_Finalize();
    return 0;
}
