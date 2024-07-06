#ifndef SYMNMFHeader
#define SYMNMFHeader

double *sym(int n, int d, double *points); 

double *ddg(int n, int d, double *points); 

double *norm(int n, int d, double *points);

double *symnmf(int n, int k, double* w, double h);

#endif