#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <pthread.h>

pthread_mutex_t mutex;

int variable = 0;


void * atStart(void * args)
{
    while (1)
    {
        printf("2 try lock\n");
        pthread_mutex_lock(&mutex);
        printf("2 locked\n");
        variable += 1;
        pthread_mutex_unlock(&mutex);
        printf("2 unlocked\n");
        sleep(1);
        printf("2 woke up\n");
    }
}


int main()
{
    pthread_t thread;
    pthread_attr_t thread_attr;
    pthread_attr_init(&thread_attr);
    pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
    pthread_mutex_init(&mutex, NULL);
    pthread_create(&thread, &thread_attr, atStart, NULL);

    sleep(3);

    printf("I START HUNTING\n");


    while (1)
    {
        printf("1 try lock\n");
        pthread_mutex_lock(&mutex);
        printf("1 locked\n");
        variable += 1;
        pthread_mutex_unlock(&mutex);
        printf("1 unlocked\n");
        usleep(500000);
        printf("%d from main", variable);
    }


    pthread_mutex_destroy(&mutex);

    return 0;
}