#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <omp.h>

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
    for (int i = 0; i < matrixHeight; ++i)
    {
        result[i] = 0;
        for (int j = 0; j < matrixWidth; ++j)
        {
            result[i] += matrix[i * matrixWidth + j] * vector[j];
        }
    }
}

void sumVectors(const double * v1, const double * v2, double * res, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        res[i] = v1[i] + v2[i];
    }
}

void subVectors(const double * v1, const double * v2, double * res, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        res[i] = v1[i] - v2[i];
    }
}

void mulVectorWithScalar(const double * vector, double * res, const int size, const double scalar)
{
    for (int i = 0; i < size; ++i)
    {
        res[i] = vector[i] * scalar;
    }
}

double countScalarSquare(const double * vector, const int size)
{
    double res = 0; 
    for (int i = 0; i < size; ++i)
    {
        res += vector[i] * vector[i];
    }
    return res;
}

int simpleIterationMethod(const double * matrix, const double * vector, double * result, const int size, const double eps, const double tao)
{
    printf("size=%d\n\n\n", size);

    double vectorBScalarSquare = countScalarSquare(vector, size);
    double * iterationVector = (double *) malloc(sizeof(double) * size); //* Ax^n - b

    double iterationVectorScalarSquare = 0;

    int iterationsCounter = 0;
    double newEps = eps * eps * vectorBScalarSquare;

    #pragma omp parallel //reduction(+ : iterationVectorScalarSquare)
    {
        int numThreads = omp_get_num_threads();
        int rankThread = omp_get_thread_num();
        int remains = size % numThreads;
        int segmentSize = size / numThreads;
        int shiftSize = segmentSize * rankThread;
        if (rankThread < remains)
        {
            segmentSize += 1;
            shiftSize += rankThread;
        }
        else
        {
            shiftSize += remains;
        }


        mulMatrixWithVector(matrix + shiftSize * size, segmentSize, size, result, iterationVector + shiftSize);

        subVectors(iterationVector + shiftSize, vector + shiftSize, iterationVector + shiftSize, segmentSize);

        double iterationVectorScalarSquarePart = countScalarSquare(iterationVector + shiftSize, segmentSize);

        #pragma omp atomic
        iterationVectorScalarSquare += iterationVectorScalarSquarePart;

        #pragma omp barrier
    
        // printf("first ivss = %lf, buf = %lf on %d\n", iterationVectorScalarSquare, buf, rankThread);

        #pragma omp barrier

        int localIterationsCounter = 0;

        while (localIterationsCounter < 10000 && iterationVectorScalarSquare >= newEps)
        {

            
            // {
            //     printf("after while x on %d with ivss = %lf\n", rankThread, iterationVectorScalarSquare);

            //     #pragma omp barrier

            //     #pragma omp master
            //     printf("========================================== %d\n", rankThread);
            // }


            // printf("%lf on iteration %d\n", countScalarSquare(iterationVector, size), iterationsCounter);
            // printf("iteration %d:\n", iterationsCounter);
            // printf("iterationVector: ");
            // printVector(iterationVector, size);
            mulVectorWithScalar(iterationVector + shiftSize, iterationVector + shiftSize, segmentSize, tao);
            // printf("iterationVector after mulVectorWithScalar: ");
            // printVector(iterationVector, size);
            // printf("result vector: ");
            // printVector(result, size);
            subVectors(result + shiftSize, iterationVector + shiftSize, result + shiftSize, segmentSize);

            #pragma omp barrier

            mulMatrixWithVector(matrix + shiftSize * size, segmentSize, size, result, iterationVector + shiftSize);
            // printf("iteration vector %d: ", iterationsCounter);
            // printVector(iterationVector, size);
            subVectors(iterationVector + shiftSize, vector + shiftSize, iterationVector + shiftSize, segmentSize);

            // printf("\n\n\nresult vector at %d iter end:\n", iterationsCounter);
            // printVector(result, size);
            // printf("\n\n\n");

            #pragma omp single
            iterationVectorScalarSquare = 0;

            iterationVectorScalarSquarePart = countScalarSquare(iterationVector + shiftSize, segmentSize);

            #pragma omp barrier

            // #pragma omp reduction(+ : iterationVectorScalarSquare)
            #pragma omp atomic
            iterationVectorScalarSquare += iterationVectorScalarSquarePart;

            // #pragma omp master
            // printf("iteration %d with ivss = %lf, that must be %lf\n", iterationsCounter, iterationVectorScalarSquare, vectorBScalarSquare * eps * eps);

            // #pragma omp barrier

            ++localIterationsCounter;

            #pragma omp barrier
        }

        #pragma omp master
        iterationsCounter = localIterationsCounter;
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



// 1678468353 = 4508

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

    //? double * buf = (double *) malloc(sizeof(double) * matrixSize);

    //? mulMatrixWithVector(A, matrixSize, matrixSize, x, buf);
    // printf("Ax = ");
    // printVector(b, matrixSize);

    // int temp = findMismatchInVectors(b, buf, matrixSize);

    // if (temp == -1) { printf("\n\nVectors are equal\n\n"); }
    // else          { printf("\n\nVectors are different in %d pos\n\n", temp); }

    //? printf("b: ");
    //? printVector(b, matrixSize);
    //? printf("Ax: ");
    //? printVector(buf, matrixSize);


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

