#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

typedef unsigned long long ull;

static const short FIRST_VECTOR_PART_MPI_ID = 0;
static const short SECOND_VECTOR_MPI_ID     = 1;
static const short ANSWER_MPI_ID            = 2;

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
    ull s = 0;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    segmentLength = (vectorLength / mpiSize) + 1;
    vectorFictiousLength = segmentLength * mpiSize;

    if (mpiRank == 0)
    {
        a = (int *) malloc(vectorFictiousLength * sizeof(int));
        b = (int *) malloc(vectorLength * sizeof(int));
        for (int i = 0; i < vectorLength; i++)
        {
            a[i] = (i * 7) % 17;
            b[i] = (vectorLength - i) * 11 % 13;
        }
        for (int i = vectorLength; i < vectorFictiousLength; i++)
        {
            a[i] = 0;
        }

        double startTime = MPI_Wtime();

        for (int i = 1; i < mpiSize; i++) 
        {
            MPI_Send(a + i * segmentLength, segmentLength, MPI_INT, i, FIRST_VECTOR_PART_MPI_ID, MPI_COMM_WORLD);
            MPI_Send(b, vectorLength, MPI_INT, i, SECOND_VECTOR_MPI_ID, MPI_COMM_WORLD);
        }

        s = countInVectors(a, b, segmentLength, vectorLength);

        ull temp;
        for (int i = 1; i < mpiSize; i++)
        {
            MPI_Recv(&temp, 1, MPI_UNSIGNED_LONG_LONG, i, ANSWER_MPI_ID, MPI_COMM_WORLD, &status);
            s += temp;
        }

        printf("%llu %lf\n", s, MPI_Wtime() - startTime);
    }
    else
    {
        a = (int *) malloc(segmentLength * sizeof(int));
        b = (int *) malloc(vectorLength  * sizeof(int));

        MPI_Recv(a, segmentLength, MPI_INT, 0, FIRST_VECTOR_PART_MPI_ID, MPI_COMM_WORLD, &status);
        MPI_Recv(b, vectorLength,  MPI_INT, 0, SECOND_VECTOR_MPI_ID,     MPI_COMM_WORLD, &status);

        s = countInVectors(a, b, segmentLength, vectorLength);

        MPI_Send(&s, 1, MPI_UNSIGNED_LONG_LONG, 0, ANSWER_MPI_ID, MPI_COMM_WORLD);
    }


    free(a);
    free(b);
    MPI_Finalize();
    return 0;
}


// 1009202029890
// 1009202029890
