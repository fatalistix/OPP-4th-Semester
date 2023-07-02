#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

typedef unsigned long long ull;

ull countInVectors(int * a, int * b, int a_size, int b_size)
{
    ull forRet = 0;
    for (int i = 0; i < a_size; i++)
    {
        for (int j = 0; j < b_size; j++)
        {
            forRet += a[i] * b[j];
        }
    }
    return forRet;
}

int main(int argc, char * argv[])
{
    int vectorLength = 145000;
    if (argc > 1)
    {
        vectorLength = atoi(argv[1]);
    }

    int * a = (int *) malloc(vectorLength * sizeof(int));
    int * b = (int *) malloc(vectorLength * sizeof(int));
    ull s = 0;
    double startTime;
    double endTime;

    for (int i = 0; i < vectorLength; i++)
    {
        a[i] = (i * 7) % 17;
        b[i] = (vectorLength - i) * 11 % 13;
    }

    MPI_Init(&argc, &argv);

    startTime = MPI_Wtime();
    s = countInVectors(a, b, vectorLength, vectorLength);
    endTime = MPI_Wtime();

    printf("%llu, %lf\n", s, endTime - startTime);

    MPI_Finalize();
    return(0);
}