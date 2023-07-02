#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>



void initFieldWithGliderTopLeftCorner(char * field, int height, int width) 
{
    memset(field, 0, width * height);
    field[0 * width + 1] = 1;
    field[1 * width + 2] = 1;
    field[2 * width + 0] = 1;
    field[2 * width + 1] = 1;
    field[2 * width + 2] = 1;
}

void printField(char * field, int height, int width)
{
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            printf("%d ", field[i * width + j]);
        }
        printf("\n");
    }
}


char countNewStageBySum(char cellState, int sum)
{
    if (cellState)
    {
        if (sum < 2 || sum > 3)
        {
            return 0;
        }
        return 1;
    }
    else
    {
        if (sum == 3)
        {
            return 1;
        }
        return 0;
    }
}

char countNewStageByNeighbours(char cellState, char one, char two, char three, char four, char five, char six, char seven, char eight)
{
    int sum = one + two + three + four + five + six + seven + eight;
    return countNewStageBySum(cellState, sum);
}

char countNewStageByIndex(char * field, int index, int fieldWidth)
{
    int sum = field[index - 1] + field[index + 1] + field[index - fieldWidth - 1] + field[index - fieldWidth] + field[index - fieldWidth + 1] + 
                    field[index + fieldWidth - 1] + field[index + fieldWidth] + field[index + fieldWidth + 1];
    return countNewStageBySum(field[index], sum);
}


void updateFieldLine(char * fieldLineSrc, char * fieldLineDst, int fieldWidth)
{
    *fieldLineDst = countNewStageByNeighbours(*fieldLineSrc,
                    *(fieldLineSrc - fieldWidth + fieldWidth - 1),
                    *(fieldLineSrc - fieldWidth),
                    *(fieldLineSrc - fieldWidth + 1),
                    *(fieldLineSrc + fieldWidth - 1),
                    *(fieldLineSrc + 1),
                    *(fieldLineSrc + fieldWidth + fieldWidth - 1),
                    *(fieldLineSrc + fieldWidth),
                    *(fieldLineSrc + fieldWidth + 1));

    for (int k = 1; k < fieldWidth - 1; ++k)
    {
        *(fieldLineDst + k) = countNewStageByIndex(fieldLineSrc, k, fieldWidth);
    }

    *(fieldLineDst + fieldWidth - 1) = countNewStageByNeighbours(*(fieldLineSrc + fieldWidth - 1),
                    *(fieldLineSrc + fieldWidth - 1 - fieldWidth - 1), 
                    *(fieldLineSrc + fieldWidth - 1 - fieldWidth), 
                    *(fieldLineSrc + fieldWidth - 1 - fieldWidth - fieldWidth + 1),
                    *(fieldLineSrc + fieldWidth - 1 - 1),
                    *(fieldLineSrc + fieldWidth - 1 - fieldWidth + 1),
                    *(fieldLineSrc + fieldWidth - 1 + fieldWidth - 1),
                    *(fieldLineSrc + fieldWidth - 1 + fieldWidth),
                    *(fieldLineSrc + fieldWidth - 1 + 1));
}




int main(int argc, char * argv[])
{
    int fieldHeight     = 100;
    int fieldWidth      = 100;
    int maxNumOfIterations = 100;

    if (argc > 1)
    {
        fieldHeight = atoi(argv[1]);
    }
    if (argc > 2)
    {
        fieldWidth  = atoi(argv[2]);
    }
    if (argc > 3)
    {
        maxNumOfIterations = atoi(argv[3]);
    }

    const int ROOT_RANK = 0;
    const int PREV_TO_NEXT_RANK_MESSAGE_ID = 101;
    const int NEXT_TO_PREV_RANK_MESSAGE_ID = 102;
    
    
    int mpiSize;
    int mpiRank;
    int mpiNextRank;
    int mpiPrevRank;

    char * field;
    int  * sendCounts;
    int  * displs;
    double start;


    char * fieldPart;
    char * fieldPartBuf;
    char * stopFlags;
    char ** previousStages;
    int currSegmentHeight;


    int minSegmentHeight;
    int heightRemains;

    MPI_Request sendToPrevReq;
    MPI_Request sendToNextReq;
    MPI_Request recvFromPrevReq;
    MPI_Request recvFromNextReq;
    MPI_Request gatherStopFlagsReq;

    MPI_Status  mpiStatus;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    minSegmentHeight = fieldHeight / mpiSize;
    heightRemains = fieldHeight - minSegmentHeight * mpiSize;
    currSegmentHeight = minSegmentHeight + ((mpiRank < heightRemains) ? 1 : 0);

    if (mpiRank == ROOT_RANK)
    {
        printf("Field height = %d\nField width  = %d\nNumber Of iterations = %d\n", fieldHeight, fieldWidth, maxNumOfIterations);
        start = MPI_Wtime();
        field = (char *) malloc(sizeof(char) * fieldWidth * fieldHeight);
        initFieldWithGliderTopLeftCorner(field, fieldHeight, fieldWidth);
        sendCounts = (int *) malloc(sizeof(int) * mpiSize);
        displs  = (int *) malloc(sizeof(int) * mpiSize);

        for (int i = 0; i < mpiSize; ++i) 
        {
            sendCounts[i] = (minSegmentHeight + ((i < heightRemains) ? 1 : 0)) * fieldWidth;
            displs[i] = (minSegmentHeight * i + ((i < heightRemains) ? i : heightRemains)) * fieldWidth;
        }
    }

    fieldPart       = (char *)  malloc(sizeof(char)   * (currSegmentHeight + 2) * fieldWidth);
    stopFlags       = (char *)  malloc(sizeof(char)   * maxNumOfIterations * mpiSize);
    previousStages  = (char **) malloc(sizeof(char *) * maxNumOfIterations);
    
    MPI_Scatterv(field, sendCounts, displs, MPI_CHAR, fieldPart + fieldWidth, currSegmentHeight * fieldWidth, MPI_CHAR, ROOT_RANK, MPI_COMM_WORLD);

    mpiPrevRank = mpiRank ? mpiRank - 1 : mpiSize - 1;
    mpiNextRank = mpiRank < mpiSize - 1 ? mpiRank + 1 : 0;
    
    int i = 0;

    for (; i < maxNumOfIterations; ++i)
    {
        MPI_Isend(fieldPart + fieldWidth, fieldWidth, MPI_CHAR, mpiPrevRank, NEXT_TO_PREV_RANK_MESSAGE_ID, MPI_COMM_WORLD, &sendToPrevReq);
        MPI_Isend(fieldPart + (currSegmentHeight) * fieldWidth, fieldWidth, MPI_CHAR, mpiNextRank, PREV_TO_NEXT_RANK_MESSAGE_ID, MPI_COMM_WORLD, &sendToNextReq);
        MPI_Irecv(fieldPart, fieldWidth, MPI_CHAR, mpiPrevRank, PREV_TO_NEXT_RANK_MESSAGE_ID, MPI_COMM_WORLD, &recvFromPrevReq);
        MPI_Irecv(fieldPart + (1 + currSegmentHeight) * fieldWidth, fieldWidth, MPI_CHAR, mpiNextRank, NEXT_TO_PREV_RANK_MESSAGE_ID, MPI_COMM_WORLD, &recvFromNextReq);
        

        for (int j = 0; j < i; ++j)
        {
            stopFlags[mpiRank * maxNumOfIterations + j] = 1;
            for (int k = 0; k < currSegmentHeight * fieldWidth; ++k)
            {
                if (previousStages[j][fieldWidth + k] != fieldPart[fieldWidth + k])
                {
                    
                    stopFlags[mpiRank * maxNumOfIterations + j] = 0;
                    break;
                }
            }
        }

        MPI_Iallgather(stopFlags + mpiRank * maxNumOfIterations, maxNumOfIterations, MPI_CHAR, stopFlags, maxNumOfIterations, MPI_CHAR, MPI_COMM_WORLD, &gatherStopFlagsReq);

        fieldPartBuf = (char *) malloc(sizeof(char) * (currSegmentHeight + 2) * fieldWidth);

        for (int k = 1; k < currSegmentHeight - 1; ++k)
        {
            updateFieldLine(fieldPart + fieldWidth * (k + 1), fieldPartBuf + fieldWidth * (k + 1), fieldWidth);
        }
        
        MPI_Wait(&sendToPrevReq, &mpiStatus);
        MPI_Wait(&recvFromPrevReq, &mpiStatus);

        updateFieldLine(fieldPart + fieldWidth, fieldPartBuf + fieldWidth, fieldWidth);
        
        MPI_Wait(&sendToNextReq, &mpiStatus);
        MPI_Wait(&recvFromNextReq, &mpiStatus);

        updateFieldLine(fieldPart + (1 + currSegmentHeight - 1) * fieldWidth, 
                        fieldPartBuf + (1 + currSegmentHeight - 1) * fieldWidth, fieldWidth);

        
        MPI_Wait(&gatherStopFlagsReq, &mpiStatus);

        // if (mpiRank == ROOT_RANK)
        // {
        //     printf("%d\n", i);
        //     printField(fieldPart + fieldWidth, currSegmentHeight, fieldWidth);
        // }

        previousStages[i] = fieldPart;
        fieldPart = fieldPartBuf;


        char finish = 0;
        for (int j = 0; j < i; ++j)
        {
            finish = 1;
            for (int k = 0; k < mpiSize; ++k)
            {
                if (!stopFlags[k * maxNumOfIterations + j])
                {
                    finish = 0;
                    break;
                }
            }
            if (finish)
            {
                break;
            }
        }


        if (finish)
        {
            break;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (mpiRank == ROOT_RANK)
    {
        printf("Took %d iterations for %lf sec for trying to return to one of previous stages\n", i + 1, MPI_Wtime() - start);
        // printField(previousStages[0], currSegmentHeight + 2, fieldWidth);
    }


    for (int j = 0; j < i; ++j)
    {
        free(previousStages[j]);
    }

    if (i != maxNumOfIterations)
    {
        free(previousStages[i]);
    }

    free(previousStages);
    free(stopFlags);
    free(fieldPart);

    if (mpiRank == ROOT_RANK)
    {
        free(field);
        free(displs);
        free(sendCounts);
    }



    MPI_Finalize();
    return 0;
}





























// fieldPartBuf[j * fieldWidth] = countNewStageByNeighbours(fieldPart[j * fieldWidth], 
//                             fieldPart[(j - 1) * fieldWidth + fieldWidth - 1], 
//                             fieldPart[(j - 1) * fieldWidth], 
//                             fieldPart[(j - 1) * fieldWidth + 1], 
//                             fieldPart[ j      * fieldWidth + fieldWidth - 1], 
//                             fieldPart[ j      * fieldWidth + 1], 
//                             fieldPart[(j + 1) * fieldWidth + fieldWidth - 1], 
//                             fieldPart[(j + 1) * fieldWidth], 
//                             fieldPart[(j + 1) * fieldWidth + 1]);
//             for (int k = 1; k < fieldWidth - 1; ++k)
//             {
//                 fieldPartBuf[j * fieldWidth + k] = countNewStageByIndex(fieldPart, j * fieldWidth + k, fieldWidth);
//             }
//             fieldPartBuf[j * fieldWidth + fieldWidth - 1] = countNewStageByNeighbours(fieldPart[j * fieldWidth + fieldWidth - 1], 
//                             fieldPart[(j - 1) * fieldWidth + fieldWidth - 1 - 1], 
//                             fieldPart[(j - 1) * fieldWidth + fieldWidth - 1], 
//                             fieldPart[(j - 1) * fieldWidth],
//                             fieldPart[ j      * fieldWidth + fieldWidth - 1 - 1],
//                             fieldPart[ j      * fieldWidth],
//                             fieldPart[(j + 1) * fieldWidth + fieldWidth - 1 - 1],
//                             fieldPart[(j + 1) * fieldWidth + fieldWidth - 1],
//                             fieldPart[(j + 1) * fieldWidth]);




// fieldPartBuf[fieldWidth] = countNewStageByNeighbours(fieldPart[fieldWidth], 
//                         fieldPart[fieldWidth - 1], 
//                         fieldPart[0], 
//                         fieldPart[1], 
//                         fieldPart[fieldWidth + fieldWidth - 1], 
//                         fieldPart[fieldWidth + 1], 
//                         fieldPart[2 * fieldWidth + fieldWidth - 1], 
//                         fieldPart[2 * fieldWidth], 
//                         fieldPart[2 * fieldWidth + 1]);
//         for (int k = 1; k < fieldWidth - 1; ++k)
//         {
//             fieldPartBuf[fieldWidth + k] = countNewStageByIndex(fieldPart, fieldWidth + k, fieldWidth);
//         }
//         fieldPartBuf[j * fieldWidth + fieldWidth - 1] = countNewStageByNeighbours(fieldPart[j * fieldWidth + fieldWidth - 1], 
//                         fieldPart[fieldWidth - 1 - 1], 
//                         fieldPart[fieldWidth - 1], 
//                         fieldPart[0],
//                         fieldPart[fieldWidth + fieldWidth - 1 - 1],
//                         fieldPart[fieldWidth],
//                         fieldPart[2 * fieldWidth + fieldWidth - 1 - 1],
//                         fieldPart[2 * fieldWidth + fieldWidth - 1],
//                         fieldPart[2 * fieldWidth]);
