#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
// C = x + A * B 
// out = x + h * dx    
double *A, *B, *C, one = 1.0;
int m,n,p; 
char *chn = "N";
A = mxGetPr(prhs[1]);
B = mxGetPr(prhs[2]);
m = mxGetM(prhs[1]);
p = mxGetN(prhs[1]);
n = mxGetN(prhs[2]);
/* Używanie jedynie w programie gdzie nie ma możliwości wprowadzenia
// złych wymiarów macierzy.
if (p != mxGetM(prhs[1])) {
mexErrMsgTxt("Inner dimensions of matrix multiply do not match");
} */
plhs[0] = mxDuplicateArray(prhs[0]);
C = mxGetPr(plhs[0]);
/* Pass all arguments by reference */
dgemm (chn, chn, &m, &n, &p, &one, A, &m, B, &p, &one, C, &m);
}