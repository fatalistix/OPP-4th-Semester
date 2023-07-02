#include "intqueue.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct iqelement_t 
{
    struct iqelement_t * next;
    int data;
} iqelement;

void iqinit(intqueue * queue) 
{
    queue->first = NULL;
    queue->last  = NULL;
    queue->size = 0;
}

int iqtop(intqueue * queue)
{
    return queue->first->data;
}

int iqempty(intqueue * queue)
{
    return !(queue->size);
}

int iqsize(intqueue * queue)
{
    return queue->size;
}

void iqpop(intqueue * queue)
{
    iqelement * firstBuf = queue->first;
    queue->first = firstBuf->next;
    free(firstBuf);
}

void iqpush(intqueue * queue, int value)
{
    if (queue->last == NULL)
    {
        queue->first = queue->last = (iqelement *) malloc(sizeof(iqelement));
    }
    else
    {
        queue->last->next = (iqelement *) malloc(sizeof(iqelement));
        queue->last = queue->last->next;
    }
    queue->last->data = value;
    queue->last->next = NULL;
}

int iqpoll(intqueue * queue)
{
    int returnData = queue->first->data;
    iqpop(queue);
    return returnData;
}

void iqfree(intqueue * queue)
{
    while (queue->first != NULL)
    {
        iqpop(queue);
    }
}