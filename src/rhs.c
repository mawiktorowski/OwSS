#include "mex.h" /* Always include this */
#include <math.h>

struct param{
    double mu, rmu, C1;
};

void rhs(double out[], double x[], double u[], struct param *p)
{
    double x1mu, x1rmu, x22, x1mu2, x1rmu2, rE, rM, rE3, rM3, coeff1, coeff2;
    
    x1mu = x[0] + p->mu;
    x1rmu = x[0] - p->rmu;
    x22 = x[1] * x[1];
    x1mu2 = x1mu * x1mu;
    x1rmu2 = x1rmu * x1rmu;
    rE = sqrt(x1mu2 + x22);
    rM = sqrt(x1rmu2 + x22);
    rE3 = rE * rE * rE;
    rM3 = rM * rM * rM;
    coeff1 = p->rmu / rE3;
    coeff2 = p->mu / rM3;
    
    out[0] = x[2];
    out[1] = x[3];
    out[2] = x[0] + 2 * x[3] - coeff1 * x1mu - coeff2 * x1rmu + u[0] * cos(u[1]) / x[4];
    out[3] = x[1] - 2 * x[2] - (coeff1 + coeff2) * x[1] + u[0] * sin(u[1]) / x[4];
    out[4] = p->C1 * u[0];
}

void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
        int nrhs, const mxArray *prhs[]) /* Input variables */
{
    struct param p;
    double *out, *x, *u;  
    mxArray *ptr;
    
    ptr = mxGetField(prhs[2], 0, "mu");
    p.mu = mxGetScalar(ptr);
    ptr = mxGetField(prhs[2], 0, "restmu");
    p.rmu = mxGetScalar(ptr);
    ptr = mxGetField(prhs[2], 0, "C1");
    p.C1 = mxGetScalar(ptr);    
    
    plhs[0] = mxCreateDoubleMatrix(1, 5, mxREAL);
 
    x = mxGetPr(prhs[0]);
    u = mxGetPr(prhs[1]);
    out = mxGetPr(plhs[0]);
    
    rhs(out, x, u, &p);
}