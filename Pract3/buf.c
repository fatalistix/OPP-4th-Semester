#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

void printVector(double * vector, int size)
{
    printf("[");
    for (int i = 0; i < size; ++i)
    {
        printf("%6.2lf, ", vector[i]);
    }
    printf("]\n");
}

double countScalarSquareParallel(double * vector, int size)
{
    double res = 0; 
    // #pragma omp for reduction(+ : res)
    for (int i = 0; i < size; ++i)
    {
        res += vector[i] * vector[i];
    }
    return res;
}

double countScalarSquare(double * vector, int size)
{
    double res = 0; 
    for (int i = 0; i < size; ++i)
    {
        res += vector[i] * vector[i];
    }
    return res;
}

void mulVectorWithScalar(const double * vector, double * res, const int size, const double scalar)
{
    for (int i = 0; i < size; ++i)
    {
        res[i] = vector[i] * scalar;
    }
}

int main()
{
    int n = 3;
    double * a = (double *) malloc(sizeof(double) * n);
    double * b = (double *) malloc(sizeof(double) * n);

    double * r = (double *) malloc(sizeof(double) * n);


    omp_set_num_threads(4);

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < 5; ++i) {
            printf("I am %d\n", omp_get_thread_num());
        }
    }

    




    return 0;
}
