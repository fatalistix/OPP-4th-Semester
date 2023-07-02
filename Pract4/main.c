/**
 * Hint for matrices' sizes:
 * 
 *           M                       K                              K
 *   |---------------|   |-------------------------|   |-------------------------|
 *   |               |   |                         |   |                         |
 *   |               |   |                         |   |                         |
 * N |               | X |         M X K           | = |                         | N
 *   |     N X M     |   |                         |   |          N X K          |
 *   |               |   |                         |   |                         |
 *   |               |   |-------------------------|   |                         |
 *   |               |                                 |                         |
 *   |---------------|                                 |-------------------------|
*/

//? 3600:3600:3600

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>


#define B_COLUMN_SEND 100
#define C_MINOR_SEND  200


void makeMatrixRandomInt(double * matrix, int height, int width, int limit)
{
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            matrix[i * width + j] = rand() % limit;
        }
    }
}

void printMatrix(double * matrix, int height, int width)
{
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            printf("%6.2lf ", matrix[i * width + j]);
        }
        printf("\n");
    }
}




int main(int argc, char * argv[])
{
    srand(12345);
    double * A = NULL;
    double * B = NULL;
    double * C = NULL;
    double * partOfA = NULL;
    

    int gridHeight = 1;
    int gridWidth  = 1;

    int N = 1;
    int M = 2;
    int K = 3;

    MPI_Comm gridComm;
    int ndims = 2;
    int dims[2];
    int periods[2] = {0, 0};
    int reorder = 1;

    MPI_Comm rowsComm;
    MPI_Comm columnsComm;
    int remainDimsRows   [2] = {0, 1};
    int remainDimsColumns[2] = {1, 0};


    int rootRank = 0;
    int rowRootRank = 0;
    int columnRootRank = 0;

    int gridRank;
    int gridCoords[2];

    double startTime;

    int segmentAHeight;
    int segmentBWidth;

    if (argc > 1)
    {
        gridHeight = atoi(argv[1]);
    }
    if (argc > 2)
    {
        gridWidth = atoi(argv[2]);
    }
    if (argc > 3)
    {
        N = atoi(argv[3]);
    }
    if (argc > 4)
    {
        M = atoi(argv[4]);
    }
    if (argc > 5)
    {
        K = atoi(argv[5]);
    }
    
    dims[0] = gridHeight;
    dims[1] = gridWidth;

    segmentAHeight = N / gridHeight;
    segmentBWidth  = M / gridWidth;

    MPI_Init(&argc, &argv);

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &gridComm);
    MPI_Cart_sub(gridComm, remainDimsRows,    &rowsComm);
    MPI_Cart_sub(gridComm, remainDimsColumns, &columnsComm);
    
    MPI_Comm_rank(gridComm, &gridRank);
    MPI_Cart_coords(gridComm, gridRank, ndims, gridCoords);


    if (gridRank == rootRank)
    {
        printf(
            "Started program with:\n"
            "grid height = %d\n"
            "grid width = %d\n"
            "first matrix size = %d\n"
            "first matrix width/second matrix height = %d\n"
            "second matrix width = %d\n", 
            gridHeight, gridWidth, N, M, K
        );

        startTime = MPI_Wtime();
        
        

        A = (double *) malloc(sizeof(double) * N * M);
        B = (double *) malloc(sizeof(double) * M * K);
        makeMatrixRandomInt(A, N, M, 10);
        makeMatrixRandomInt(B, M, K, 20);


        // printMatrix(A, N, M);
        // printf("==============================================\n");
        // printMatrix(B, M, K);
        // printf("==============================================\n");

        MPI_Datatype sendColumn;
        MPI_Type_vector(M, segmentBWidth, K, MPI_DOUBLE, &sendColumn);
        MPI_Type_commit(&sendColumn);

        for (int i = 1; i < gridWidth; ++i)
        {
            MPI_Send(B + i * segmentBWidth, 1, sendColumn, i, B_COLUMN_SEND + i, gridComm);
        }

        MPI_Bcast(B, 1, sendColumn, columnRootRank, columnsComm);

        MPI_Type_free(&sendColumn);
    }
    else
    {
        MPI_Datatype recvContiguos;
        MPI_Type_contiguous(M * segmentBWidth, MPI_DOUBLE, &recvContiguos);
        MPI_Type_commit(&recvContiguos);
        B = (double *) malloc(sizeof(double) * M * segmentBWidth);
        MPI_Status status;
        if (gridCoords[0] == 0)
        {   
            MPI_Recv(B, 1, recvContiguos, rootRank, B_COLUMN_SEND + gridRank, gridComm, &status);
        }

        MPI_Bcast(B, 1, recvContiguos, columnRootRank, columnsComm);


        MPI_Type_free(&recvContiguos);
    }


    


    partOfA = (double *) malloc(sizeof(double) * segmentAHeight * M);


    if (gridCoords[1] == 0)
    {
        MPI_Scatter(A, segmentAHeight * M, MPI_DOUBLE, partOfA, segmentAHeight * M, MPI_DOUBLE, columnRootRank, columnsComm);
    }



    MPI_Bcast(partOfA, segmentAHeight * M, MPI_DOUBLE, rowRootRank, rowsComm);





    if (gridRank == rootRank)
    {
        C = (double *) calloc(sizeof(double), N * K);
    }
    else 
    {
        C = (double *) calloc(sizeof(double), segmentAHeight * segmentBWidth);
    }



    if (gridRank == rootRank)
    {
        for (int i = 0; i < segmentAHeight; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                for (int k = 0; k < segmentBWidth; ++k)
                {
                    C[i * K + k] += partOfA[i * M + j] * B[j * K + k];
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < segmentAHeight; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                for (int k = 0; k < segmentBWidth; ++k)
                {
                    C[i * segmentBWidth + k] += partOfA[i * M + j] * B[j * segmentBWidth + k];
                }
            }
        }
    }


    

    
    if (gridRank == rootRank)
    {
        MPI_Datatype minor;
        MPI_Type_vector(segmentAHeight, segmentBWidth, K, MPI_DOUBLE, &minor);
        MPI_Type_commit(&minor);
        MPI_Status status;
        for (int i = 0; i < gridHeight; ++i)
        {
            for (int j = 0; j < gridWidth; ++j)
            {
                if (i != 0 || j != 0)
                {
                    MPI_Recv(C + i * segmentAHeight * K + j * segmentBWidth, 1, minor, i * gridWidth + j, C_MINOR_SEND + i * gridWidth + j, gridComm, &status);
                }
            }
        }
        MPI_Type_free(&minor);

        // printMatrix(C, N, K);
    }
    else
    {
        MPI_Send(C, segmentAHeight * segmentBWidth, MPI_DOUBLE, rootRank, C_MINOR_SEND + gridRank, gridComm);
        // printf("%d\n", gridRank);
    }

    
    if (gridRank == rootRank) {
        printf("Answer found for %lf secs\n", MPI_Wtime() - startTime);
        free(A);
    } 
    free(partOfA);
    free(B);
    free(C);




    MPI_Comm_free(&gridComm);
    MPI_Comm_free(&columnsComm);
    MPI_Comm_free(&rowsComm);

    MPI_Finalize();
    return 0;
}










































































































































// int main(int argc, char * argv[]) 
// {
//     int * A;
//     int * partOfA;
//     int rank;
//     MPI_Init(&argc, &argv);
//     partOfA = (int *) malloc(sizeof(int) * 6 * 6);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Datatype col, coltype;
//     MPI_Type_vector(6, 2, 7, MPI_INT, &col);
//     MPI_Type_commit(&col);

//     MPI_Datatype col_resized;
//     MPI_Type_create_resized(col, 0, 6, &col_resized);
//     MPI_Type_commit(&col_resized);



//     if (rank == 0)
//     {
//         A = (int *) malloc(sizeof(int) * 6 * 6);
//         for (int i = 0; i < 36; ++i) { A[i] = i; }
//         MPI_Send(A + 1, 1, col_resized, 1, 0, MPI_COMM_WORLD);
//     }

//     if (rank == 1)
//     {
//         MPI_Status stat;
//         for (int i = 0; i < 36; ++i) { partOfA[i] = 0; }
//         MPI_Recv(partOfA, 1, col_resized, 0, 0, MPI_COMM_WORLD, &stat);
//         for (int i = 0; i < 6; ++i)
//         {
//             for (int j = 0; j < 6; ++j)
//             {
//                 printf("%d ", partOfA[i * 6 + j]);
//             }
//             printf("\n");
//         }
//     }



//     MPI_Finalize();
//     return 0;
// }






























































// int main(int argc, char * argv[]) 
// {
//     int hGridSize;
//     int vGridSize;

//     if (argc < 3)
//     {
//         fprintf(stderr, "Expected 2 arguments: grid's height and grid's width");
//         return EXIT_FAILURE;
//     }

//     hGridSize = atoi(argv[1]);
//     vGridSize = atoi(argv[2]);

//     MPI_Init(&argc, &argv);
//     MPI_Comm grid_comm;
//     int dims[2] = {hGridSize, vGridSize};
//     int periods[2] = {0, 0};

//     MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);


//     int coords[2];
//     int grid_comm_rank;
//     MPI_Comm_rank(grid_comm, &grid_comm_rank);
//     MPI_Cart_coords(grid_comm, grid_comm_rank, 2, coords);

//     int grid_comm_size;
//     MPI_Comm_size(grid_comm, &grid_comm_size);

//     // printf("%d %d on rank %d of %d\n", coords[0], coords[1], grid_comm_rank, grid_comm_size);
//     coords[1]++;
//     MPI_Cart_rank(grid_comm, coords, &grid_comm_rank);
//     printf("%d rank on %d %d\n", grid_comm_rank, coords[0], coords[1]);

//     MPI_Finalize();
//     return 0;
// }

