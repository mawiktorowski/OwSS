#include "mex.h" /* Always include this */
#include <math.h>

struct param{
    double gE, gM, D, omega, C1;
};

void rhs(double out[], double t, double x[], double u[], struct param *p)
{
    double xm, ym, xm2, ym2, x12, x22, rE, rM, rE3, rM3, coeff1, coeff2; 
    
    xm = x[0] - p->D * cos(p->omega * t);
    ym = x[1] - p->D * sin(p->omega * t);
    xm2 = xm * xm;
    ym2 = ym * ym;
    x12 = x[0] * x[0];
    x22 = x[1] * x[1];
    rE = sqrt(x12 + x22);
    rM = sqrt(xm2 + ym2);
    rE3 = rE * rE * rE;
    rM3 = rM * rM * rM; 
    coeff1 = p->gE / rE3;
    coeff2 = p->gM / rM3;
    
    out[0] = x[2];
    out[1] = x[3];
    out[2] = - coeff1 * x[0] - coeff2 * xm + u[0] * cos(u[1]) / x[4];
    out[3] = - coeff1 * x[1] - coeff2 * ym + u[0] * sin(u[1]) / x[4];
    out[4] = p->C1 * u[0];
}

void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
        int nrhs, const mxArray *prhs[]) /* Input variables */
{
    struct param p;
    double *out, t, *x, *u;  
    
    p.gE = mxGetScalar(mxGetField(prhs[3], 0, "gE"));
    p.gM = mxGetScalar(mxGetField(prhs[3], 0, "gM"));
    p.D = mxGetScalar(mxGetField(prhs[3], 0, "D"));
    p.omega = mxGetScalar(mxGetField(prhs[3], 0, "omega"));    
    p.C1 = mxGetScalar(mxGetField(prhs[3], 0, "C1"));  
    
    plhs[0] = mxCreateDoubleMatrix(1, 5, mxREAL);
 
    t = mxGetScalar(prhs[0]);
    x = mxGetPr(prhs[1]);
    u = mxGetPr(prhs[2]);
    out = mxGetPr(plhs[0]);
    
    rhs(out, t, x, u, &p);
}