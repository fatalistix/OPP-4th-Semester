#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pthread.h>

static const int ANY_SOURCE_LISTEN = 100;
static const int SHARE_SIZE = 102;
static const int SHARE_TASKS = 103;


static const int MESSAGE_NO_TASKS = -101;


static const double SHARE_BYPASS = 1./10.;
static const double SHARE_PART   = 2./3.;


static const int STOP_LISTEN = -1;



int * tasksInputData;
int taskDataSize = 500;
int numOfTasks;
double globalRes = 0;
int iterCounter = 0;
int numOfCompleted;

pthread_mutex_t mutex;





int completeTasks()
{
    while (numOfCompleted < numOfTasks)
    {
        pthread_mutex_lock(&mutex);
        int repeatNum = tasksInputData[numOfCompleted];
        ++numOfCompleted;
        pthread_mutex_unlock(&mutex);
        for (int j = 0; j < repeatNum; ++j)
        {
            globalRes += sin(j);
        }
    }
    return numOfCompleted;
}


char canSend()
{
    return (numOfTasks - numOfCompleted) > taskDataSize * SHARE_BYPASS;
}


void * listener(void * args)
{
    int askerRank;
    MPI_Status status;
    while (1)
    {
        MPI_Recv(&askerRank, 1, MPI_INT, MPI_ANY_SOURCE, ANY_SOURCE_LISTEN, MPI_COMM_WORLD, &status);
        if (askerRank != STOP_LISTEN)
        {
            pthread_mutex_lock(&mutex);
            if (canSend())
            {
                int remains = numOfTasks - numOfCompleted;
                int sendArraySize = remains - (int) (remains * SHARE_PART);
                numOfTasks -= sendArraySize;
                pthread_mutex_unlock(&mutex);

                 
                MPI_Send(&sendArraySize, 1, MPI_INT, askerRank, SHARE_SIZE, MPI_COMM_WORLD);
                MPI_Send(tasksInputData + numOfCompleted, sendArraySize, MPI_INT, askerRank, SHARE_TASKS, MPI_COMM_WORLD);
               
            }
            else 
            {
                MPI_Send(&MESSAGE_NO_TASKS, 1, MPI_INT, askerRank, SHARE_SIZE, MPI_COMM_WORLD);
                pthread_mutex_unlock(&mutex);
            }
        }
        else 
        {
            break;
        }
    }

    pthread_exit(NULL);
}


int main(int argc, char * argv[])
{
    int iterCounterMax = 3;
    int taskFactor = 10000;


    if (argc > 1)
    {
        iterCounterMax = atoi(argv[1]);
    }
    if (argc > 2)
    {
        taskFactor = atoi(argv[2]);
    }
    if (argc > 3)
    {
        taskDataSize = atoi(argv[3]);
    }

    int mpiSize;
    int mpiRank;

    int requiedSupport = MPI_THREAD_MULTIPLE;
    int realSupport;

    MPI_Init_thread(&argc, &argv, requiedSupport, &realSupport);
    if (realSupport != MPI_THREAD_MULTIPLE)
    {
        fprintf(stderr, "%s: MPI_THREAD_MULTIPLE is not supported\n", argv[0]);
        return 1;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if (mpiRank == 0)
    {
        printf("%s: running with iterCounterMax = %d, taskFactor = %d, taskDataSize = %d\n", argv[0], iterCounterMax, taskFactor, taskDataSize);
        printf("<===================================================>\n");
    }

    tasksInputData = (int *) malloc(sizeof(int) * taskDataSize);

    
    pthread_t listenThread;
    pthread_attr_t attr;
    if (pthread_attr_init(&attr))
    {
        perror("Error loading attributes");
        MPI_Abort(MPI_COMM_WORLD, errno);
    }
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE))
    {
        perror("Cannot apply JOINABLE for attr");
        MPI_Abort(MPI_COMM_WORLD, errno);
    }

    if (pthread_mutex_init(&mutex, NULL))
    {
        perror("Cannot init mutex");
        MPI_Abort(MPI_COMM_WORLD, errno);
    }

    if (pthread_create(&listenThread, &attr, listener, NULL))
    {   
        perror("Cannot create listener thread");
        MPI_Abort(MPI_COMM_WORLD, errno);
    }

    double startCountTime;
    double currRankCountTime;
    double maxCountTime;
    double minCountTime;
    double maxTimeDiff;
    MPI_Status status;

    for (; iterCounter < iterCounterMax; ++iterCounter)
    {
        int tasksCompleted = 0;
        numOfTasks = taskDataSize;
        startCountTime = MPI_Wtime(); 
        for (int i = 0; i < taskDataSize; ++i)
        {
            tasksInputData[i] = abs(taskDataSize / 2 - (i + taskDataSize * mpiRank) % 100) * abs(mpiRank - (iterCounter % mpiSize)) * taskFactor;
        }
        
        tasksCompleted += completeTasks();

        printf("%d completed %d\n", mpiRank, tasksCompleted);

        for (int i = 0; i < mpiSize; ++i)
        {
            if (mpiRank == i)
            {
                continue;
            }
            int otherTaskArraySize;

            MPI_Send(&mpiRank, 1, MPI_INT, i, ANY_SOURCE_LISTEN, MPI_COMM_WORLD);
            MPI_Recv(&otherTaskArraySize, 1, MPI_INT, i, SHARE_SIZE, MPI_COMM_WORLD, &status);
            printf("%d can took %d from %d\n", mpiRank, otherTaskArraySize, i);
            if (otherTaskArraySize == MESSAGE_NO_TASKS)
            {
                continue;
            }
            MPI_Recv(tasksInputData, otherTaskArraySize, MPI_INT, i, SHARE_TASKS, MPI_COMM_WORLD, &status);
            pthread_mutex_lock(&mutex);
            numOfCompleted = 0;
            numOfTasks = otherTaskArraySize;
            pthread_mutex_unlock(&mutex);
            tasksCompleted += completeTasks();
        }

        currRankCountTime = MPI_Wtime() - startCountTime;

        MPI_Allreduce(&currRankCountTime, &maxCountTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&currRankCountTime, &minCountTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);


        maxTimeDiff = maxCountTime - minCountTime;


        printf("%d: %d task completed during iteration: %d; result: %lf; TOOK: %lf; DELTA: %lf; UNBALANCE: %lfPERC\n", mpiRank, tasksCompleted, iterCounter, globalRes, currRankCountTime, maxTimeDiff, maxTimeDiff / maxCountTime * 100);
    }

    MPI_Send(&STOP_LISTEN, 1, MPI_INT, mpiRank, ANY_SOURCE_LISTEN, MPI_COMM_WORLD);



    if (pthread_attr_destroy(&attr))
    {
        perror("Error destroying attr");
    }

    if (pthread_mutex_destroy(&mutex))
    {
        perror("Cannot destroy mutes");
    }

    if (pthread_join(listenThread, NULL))
    {
        perror("Error joining thread");
    }

    free(tasksInputData);

    MPI_Finalize();
    return 0;
}