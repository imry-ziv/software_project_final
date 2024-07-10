#ifndef SYMNMFHeader
#define SYMNMFHeader

double *sym(int n, int d, double *points, int *status); 

double *ddg(int n, int d, double *points, int *status); 

double *norm(int n, int d, double *points, int *status);

double *symnmf(int n, int k, double* w, double *h, int *status);

#endif