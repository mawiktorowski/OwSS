#include "mex.h" /* Always include this */
#include <math.h>

struct param{
    double mu, rmu, C1;
};

void rhsPsi(double out[], double x[], double u[], struct param *p)
{
    
    double x1mu, x1rmu, x22, x1mu2, x1rmu2, rE, rM, rE3, rM3, coeff1, coeff2;
    double rE2, rM2, d31, d41, d32, d42, x52, sphi, cphi, coeff3, coeff4;
    
    x1mu = x[0] + p->mu;
    x1rmu = x[0] - p->rmu;
    x22 = x[1] * x[1];
    x1mu2 = x1mu * x1mu;
    x1rmu2 = x1rmu * x1rmu;
    rE = sqrt(x1mu2 + x22);
    rM = sqrt(x1rmu2 + x22);
    rE2 = rE * rE;
    rM2 = rM * rM;
    rE3 = rE2 * rE;
    rM3 = rM2 * rM;    
    coeff1 = p->rmu / rE3;
    coeff2 = p->mu / rM3;    
    coeff3 = coeff1 / rE2;
    coeff4 = coeff2 / rM2;
    
    x52 = x[4] * x[4];
    sphi = sin(u[1]);
    cphi = cos(u[1]);
    
    d31 = - 1 + coeff3 * (rE2 - 3 * x1mu2) + coeff4 * (rM2 - 3 * x1rmu2);
    d41 = - 3 * x[1] * (coeff3 * x1mu + coeff4 * x1rmu);
    d32 = d41; 
    d42 = - 1 + coeff3 * (rE2 - 3 * x22) + coeff4 * (rM2 - 3 * x22);
    
    out[0] = x[2];
    out[1] = x[3];
    out[2] = x[0] + 2 * x[3] - coeff1 * x1mu - coeff2 * x1rmu + u[0] * cphi / x[4];
    out[3] = x[1] - 2 * x[2] - (coeff1 + coeff2) * x[1] + u[0] * sphi / x[4];    
    out[4] = p->C1 * u[0];
    // rownania sprzezone
    out[5] = x[7] * d31 + x[8] * d41;    
    out[6] = x[7] * d32 + x[8] * d42;    
    out[7] = - x[5] + 2 * x[8];   
    out[8] = - x[6] - 2 * x[7]; 
    out[9] = u[0] * (x[7] * cphi + x[8] * sphi) / x52;
    // dH/du
    out[10] = x[7] * cphi / x[4] + x[8] * sphi / x[4] + x[9] * p->C1;          
    out[11] = u[0] * (- x[7] * sphi + x[8] * cphi) / x[4];        
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
    
    plhs[0] = mxCreateDoubleMatrix(1, 12, mxREAL);
    
    x = mxGetPr(prhs[0]);
    u = mxGetPr(prhs[1]);
    out = mxGetPr(plhs[0]);
    
    rhsPsi(out, x, u, &p);
}
