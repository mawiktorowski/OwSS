#include "mex.h" /* Always include this */
#include <math.h>

struct param{
    double gE, gM, D3, C1;
};

void rhs(double out[], double x[], double u[], struct param *p)
{
    double xm, ym, xm2, ym2, x32, x42, rE, rM, rE3, rM3,
            coeff1, coeff2, coeff3; 
    
    xm = x[2] - x[0];
    ym = x[3] - x[1];
    xm2 = xm * xm;
    ym2 = ym * ym;
    x32 = x[2] * x[2];
    x42 = x[3] * x[3];
    rE = sqrt(x32 + x42);
    rM = sqrt(xm2 + ym2);
    rE3 = rE * rE * rE;
    rM3 = rM * rM * rM;
    coeff1 = p->gE / p->D3;    
    coeff2 = p->gE / rE3;
    coeff3 = p->gM / rM3;
    
    out[0] = x[4];
    out[1] = x[5];
    out[2] = x[6];
    out[3] = x[7];
    out[4] = - coeff1 * x[0];
    out[5] = - coeff1 * x[1];
    out[6] = - coeff2 * x[2] - coeff3 * xm + u[0] * cos(u[1]) / x[8];
    out[7] = - coeff2 * x[3] - coeff3 * ym + u[0] * sin(u[1]) / x[8];
    out[8] = p->C1 * u[0];
}

void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
        int nrhs, const mxArray *prhs[]) /* Input variables */
{
    struct param p;
    double *out, *x, *u;  
    
    p.gE = mxGetScalar(mxGetField(prhs[2], 0, "gE"));
    p.gM = mxGetScalar(mxGetField(prhs[2], 0, "gM"));
    p.D3 = mxGetScalar(mxGetField(prhs[2], 0, "D3"));  
    p.C1 = mxGetScalar(mxGetField(prhs[2], 0, "C1"));  
    
    plhs[0] = mxCreateDoubleMatrix(1, 9, mxREAL);
 
    x = mxGetPr(prhs[0]);
    u = mxGetPr(prhs[1]);
    out = mxGetPr(plhs[0]);
    
    rhs(out, x, u, &p);
}