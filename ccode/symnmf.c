#ifndef SYMNMFSource
#define SYMNMFSource

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "symnmf.h"
#include <string.h>

#define beta 0.5
#define eps 1e-4
#define maxIteration 300
#define PRINTERROR printf("An Error Has Occured\n")


/*
Implementation convention - matrices are stored as linear arrays.
Therefore points [d,n] will be an d*n array where [0], ..., [d-1] 
comprise the first vector, etc. 
Regarding memory - method do not clean input memory, and 
only clean allocated memory that's not the result.
*/
int Index(int rowLength, int row, int col)
{
    return row * rowLength + col;
}

double EuclideanDistanceSqr(int d, double *points, int i1, int i2)
{
    double sum = 0, temp = 0;
    int i;
    for (i = 0; i < d; i++)
    {
        temp = points[i1+i] - points[i2+i];
        sum += temp * temp;
    }
    return sum;
}

double *AllocateMatrix(int a, int b, int* status)
{
    double *res = (double*)malloc(a * b * sizeof(double));
    if (NULL == res)
        *status = 1;
    return res;
}

/*
A[a,b], B[b,c], and dest[a,c]
overwrites result into dest
*/
void MultiplyMatricesNonAlloc(int a, int b, int c, double *A, double *B, double *dest)
{
    double sum;
    int i, j, k;
    for (i = 0; i < a; i++)
    {
        for (j = 0; j < c; j++)
        {
            sum = 0;
            for (k = 0; k < b; k++)
            {
                sum += A[Index(a, i, k)] * B[Index(b, k, j)];
            }
            dest[Index(a, i, j)] = sum;
        }
    }
}

/* A[a,b] and B[b,c] */
double *MultiplyMatrices(int a, int b, int c, double *A, double *B, int* status)
{
    double *res = AllocateMatrix(a, c, status);
    if (0 != *status)
        return NULL;
    MultiplyMatricesNonAlloc(a, b, c, A, B, res);
    return res;
}

double *sym(int n, int d, double *points, int* status)
{
    double dist;
    int i, j;
    double *A = AllocateMatrix(n, n, status);
    printf("Entering sym\n");
    if (0 != *status)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
                A[Index(n, i, j)] = 0.0;
            else 
            {
                dist = EuclideanDistanceSqr(d, points, Index(d, i, 0), Index(d, j, 0));
                A[Index(n, i, j)] = exp(-dist / 2.0);
            }
        }
    }
    printf("Exiting sym\n");
    return A;
}

double *ddgFromA(int n, double *A, int* status)
{    
    int i, j;
    double sum;
    double *D = AllocateMatrix(n, n, status);
    if (0 != *status)
        return NULL;
    for (i = 0; i < n; i++)
    {
        sum = 0;
        for (j = 0; j < n; j++)
        {
            sum += A[Index(n, i, j)];
            D[Index(n, i, j)] = 0.0;
        }
        D[Index(n, i, i)] = sum;
    }
    return D;
}

double *ddg(int n, int d, double *points, int* status)
{
    double *res, *D = sym(n, d, points, status);
    if (0 != *status)
    {
        printf("Failed ddg\n");
        return NULL;
    }
    res = ddgFromA(n, D, status);
    if (0 != *status)
    {
        printf("ddg failed in ddgFromA\n");
        PRINTERROR;
    }
    free(D);
    return res; /* can be null if we failed along the way */
}

double *norm(int n, int d, double *points, int* status)
{
    int i;
    double *A, *D, *temp;
    A = sym(n, d, points, status);
    if (0 != *status)
        return NULL;
    D = ddgFromA(n, A, status);
    if (0 != *status)
        goto norm_free1;
    for (i = 0; i < n; i++)
    {
        D[Index(n,i,i)] = 1.0/sqrt(D[Index(n,i,i)]);
    }
    temp = MultiplyMatrices(n, n, n, D, A, status);
    if (0 != *status)
        goto norm_free2;
    /*
    A is overwritten here since it's no longer in use
    It now stores the result of the calculation.
    */
    MultiplyMatricesNonAlloc(n, n, n, temp, D, A);
norm_free2:
    free(temp); 
norm_free1:
    free(D);
    return A;
}

float F2NormOfDifference(int a, int b, double *A, double *B)
{
    int i, j, index;
    double dif, sum = 0;
    for (i = 0; i < a; i++)
    {
        for (j = 0; j < b; j++)
        {
            index = Index(a, i, j);
            dif = A[index] - B[index];
            sum += dif * dif;
        }
    }
    return sum;
}

/* src[a,b], dest[b,a] */
void Transpose(int a, int b, double *src, double *dest)
{
    int i, j;
    for (i = 0; i < a; i++)
    {
        for (j = 0; j < b; j++)
        {
            dest[Index(a, j, i)] = src[Index(b, i, j)];
        }
    }
}

/* H[n,k] (cur and next), W[n,n], temp1[n*k], temp2[k,k]
since we use temp1 as both [n,k]  and [k,n]
*/
void ConvergenceStep(int n, int k, double* W, double *Hcur, double *Hnext, double *temp1, double *temp2)
{
    int i,j, index;
    /* Step 1: calculate W*Hcur
    we'll store this in Hnext since it's [n,k] */
    MultiplyMatricesNonAlloc(n, n, k, W, Hcur, Hnext);
    /* Step 2: calculate Hcur^T * Hcur
    We'll store H^T in temp1, and the resulting [k,k] matrix in temp2 */
    Transpose(n, k, Hcur, temp1);
    MultiplyMatricesNonAlloc(k, n, k, temp1, Hcur, temp2);
    /*Step 3: calculate Hcur * temp2
    We'll store the result back into temp1 */
    MultiplyMatricesNonAlloc(n, k, k, Hcur, temp2, temp1);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            index = Index(k, i, j);
            Hnext[index] = Hcur[index] * (1 - beta + beta * Hnext[index] / temp1[index]);
        }
    }
}

double *symnmf(int n, int k, double* W, double *H, int* status)
{
    double *Hcur, *Hnext, *temp1, *temp2, *swap;
    int i;
    Hcur = H;
    Hnext = AllocateMatrix(n, k, status);
    if (0 != *status)
        goto symnmf_free1;
    temp1 = AllocateMatrix(n, k, status);
    if (0 != *status)
        goto symnmf_free2;
    temp2 = AllocateMatrix(k, k, status);
    if (0 != *status)
        goto symnmf_free3;
    for (i = 0; i < maxIteration; i++)
    {
        ConvergenceStep(n, k, W, Hcur, Hnext, temp1, temp2);
        swap = Hcur;
        Hcur = Hnext;
        Hnext = swap;
        if (F2NormOfDifference(n, k, Hcur, Hnext) < eps)
            break;
    }
symnmf_free3:
    free(temp2); 
symnmf_free2:
    free(temp1); 
symnmf_free1:
    free(Hnext); 
    return Hcur;
}


void PrintPoint(int d, int index, double *points)
{
    int elem;
    printf("%.4f", points[index]);
    for (elem = 1; elem < d; elem++)
    {
        printf(",%.4f", points[index + elem]);
    }
    printf("\n");
}

void PrintPoints(int n, double *points)
{
	int i;
	for (i = 0; i < n; i++)
	{
		PrintPoint(n, n * i, points);
	}
}

void ReadPoints(FILE *stream, double *points, int entryCount, int *d, int *n, int* status)
{
    int convs, pointIndex, elem;
    char sep;
    elem = 0;
    while (1)
    {
        convs = fscanf(stream, "%lf", points + elem);
        ++elem;
        sep = fgetc(stream);
        if (sep == ',')
            continue;
        if (sep == '\n')
            break;
        *status = 1;
        return; 
    }
    *d = elem;
    *n = entryCount / *d;
    printf("Dimension: %d, number of points: %d", *d, *n);
    if ((*n) * (*d) != entryCount)
    {
        *status = 1;
        PRINTERROR;
        return;
    }
    pointIndex = 1;
    for (pointIndex = 1; pointIndex < *n; ++pointIndex)
	{
		for (elem = 0; elem < *d; elem++)
		{
			convs = fscanf(stream, "%lf", points + Index(*d, pointIndex, elem));
			if (1 != convs)
            {
				*status = 1;
                PRINTERROR;
                return;
            }
			sep = fgetc(stream);
			if (sep == ',' && elem < *d - 1)
				continue;
			if (sep == '\n' && elem == *d - 1)
				continue;
            *status = 1;
            PRINTERROR;
            return;
		}
	}
    /* validate EOF - they say we dont actually need to do this
    sep = fgetc(stream);
    if (sep != EOF)
    {
        *status = 1;
        printf("Did not reach EOF: instead got %d after %d points.\n", sep, pointIndex);
        sep = fgetc(stream);
        while (sep != EOF)
        {
            printf("%c\n", sep = fgetc(stream));
        } 
        PRINTERROR;
        return;
    }
    */
}

int CountPoints(FILE *stream)
{/* each point has a decimal point, so we'll count those */
    int res;
    char c;
    res = 0;
    do
    {
        c = fgetc(stream);
        if (c == '.')
            ++res;
    } while (c != EOF);
    rewind(stream);
    return res;
}

int main(int argc, char **argv)
{
    int n, d, status, entryCount; 
    char *goal;
    double *points, *res;
    FILE *stream;
    if (3 != argc)
    {
        PRINTERROR;
        return 1;
    }
    goal = argv[1];
    stream = fopen(argv[2], "r");
    status = 0;
    /* dont let your memes be dreams */
    points = AllocateMatrix(entryCount = CountPoints(stream), 1, &status);
    if (0 != status)
        return 1;
    n = 0;
    d = 0;
    ReadPoints(stream, points, entryCount, &d, &n, &status);
    fclose(stream);
    if (0 != status)
        goto main_free1;
    if (0 == strcmp(goal, "sym"))
    {
        res = sym(n, d, points, &status);
    }
    else if (0 == strcmp(goal, "ddg"))
    {
        printf("Enterind ddg\n");
        res = ddg(n, d, points, &status);
    }
    else if (0 == strcmp(goal, "norm"))
    {
        res = norm(n, d, points, &status);
    }
    if (NULL != res)
    {
        PrintPoints(n, res);
        free(res);
    }
main_free1:
    free(points);
    return status;
}

#endif