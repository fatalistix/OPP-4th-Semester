#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char * argv[]) 
{
    int * A;
    int * partOfA;
    int rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // здесь тип вектор
    MPI_Datatype col, coltype;
    MPI_Type_vector(6, 2, 6, MPI_INT, &col);
    MPI_Type_commit(&col);

    // MPI_Datatype col_resized;
    // MPI_Type_create_resized(col, 0, 6, &col_resized);
    // MPI_Type_commit(&col_resized);

    partOfA = (int *) malloc(sizeof(int) * 6 * 6);

    

    if (rank == 0)
    {
        // здесь заполняется матрица как сверху слева в терминале
        A = (int *) malloc(sizeof(int) * 6 * 6);
        for (int i = 0; i < 36; ++i) { A[i] = i; }


        // отправляю свой тип данных в единичном количестве первому процессу (0 => 1)
        MPI_Send(A + 1, 1, col, 1, 0, MPI_COMM_WORLD);


        // вывод А
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                printf("%3d ", A[i * 6 + j]);
            }
            printf("\n");
        }
        printf("==========================================\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 1)
    {
        MPI_Status stat;
        // зануляю матрицу перед приемом
        for (int i = 0; i < 36; ++i) { partOfA[i] = 0; }


        // мкняю значения получения на набор даблов
        MPI_Recv(partOfA, 12, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);

        // вывод результата
        // они пришли подряд, хотя я отправлял вектор, а приняо дабл
        
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                printf("%3d ", partOfA[i * 6 + j]);
            }
            printf("\n");
        }
    }



    MPI_Finalize();
    return 0;
}