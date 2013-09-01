#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    double *A, *B, *C;
    int m, n, i;
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
// Brak sprawdzania błędów
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    C = mxGetPr(plhs[0]);
    /* Pass all arguments by reference */
    for (i = 0; i < m * n; i++)
    {
        C[i] = A[i] + B[i];
    }
}
