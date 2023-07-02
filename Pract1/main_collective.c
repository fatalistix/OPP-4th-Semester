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
    int vectorLength = argc > 1 ? atoi(argv[1]) : 145000;
    int vectorFictiousLength;
    int mpiSize;
    int mpiRank;
    int segmentLength;
    int * a;
    int * b;
    int * aPartCopy;
    ull localAnswer = 0;
    ull answer = 0;

    MPI_Status status;

    double startTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    segmentLength = (vectorLength / mpiSize) + 1;
    vectorFictiousLength = segmentLength * mpiSize;

    b = (int *) malloc(vectorLength  * sizeof(int));

    if (mpiRank == 0)
    {
        a = (int *) malloc(vectorFictiousLength * sizeof(int));
        for (int i = 0; i < vectorLength; i++)
        {
            a[i] = (i * 7) % 17;
            b[i] = ((vectorLength - i) * 11) % 13;
        }
        for (int i = vectorLength; i < vectorFictiousLength; i++)
        {
            a[i] = 0;
        }
        startTime = MPI_Wtime();
    }
    
    aPartCopy = (int *) malloc(segmentLength * sizeof(int));

    MPI_Scatter(a, segmentLength, MPI_INT, aPartCopy, segmentLength, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast  (b, vectorLength,  MPI_INT, 0, MPI_COMM_WORLD);

    localAnswer = countInVectors(aPartCopy, b, segmentLength, vectorLength);
    
    MPI_Reduce(&localAnswer, &answer, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);


    if (mpiRank == 0)
    {
        printf("%llu %lf\n", answer, MPI_Wtime() - startTime);
        free(a);
    }

    free(aPartCopy);
    free(b);
    MPI_Finalize();
    return 0;
}