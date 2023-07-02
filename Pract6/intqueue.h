#pragma once

typedef struct iqelement_t iqelement;

typedef struct intqueue_t 
{
    iqelement * first;
    iqelement * last;
    int size;
} intqueue;

void iqinit (intqueue * queue);
int  iqtop  (intqueue * queue);
int  iqempty(intqueue * queue);
int  iqsize (intqueue * queue);
void iqpop  (intqueue * queue);
void iqpush (intqueue * queue, int value);
int  iqpoll (intqueue * queue);
void iqfree (intqueue * queue);
