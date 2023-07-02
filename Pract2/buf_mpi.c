#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

int main(int argc, char * argv[])
{
    int * arr = (int *) malloc(sizeof(int) * 9);
    int * res = (int *) calloc(sizeof(int), 9);
    int * par = (int *) malloc(sizeof(int) * 3);
    int * pos = (int *) malloc(sizeof(int) * 3);

    int rank;

    for (int i = 0; i < 9; ++i)
    {
        arr[i] = i + 10;
    }

    par[0] = 4;
    par[1] = 3;
    par[2] = 2;

    pos[0] = 0;
    pos[1] = 4;
    pos[2] = 4 + 3;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Scatterv(arr, par, pos, MPI_INT, res, 4, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 2)
    {
        for (int i = 0; i < 9; ++i)
        {
            printf("%d ", res[i]);
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}