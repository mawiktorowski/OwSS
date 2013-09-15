#include "mex.h" /* Always include this */
#include <math.h>

struct param{
    double gE, gM, D, omega, C1;
};

void rhsPsi(double out[], double t ,double x[], double u[], struct param *p)
{
    double xm, ym, xm2, ym2, x12, x22, rE, rM, rE2, rM2, rE3, rM3,
            coeff1, coeff2, coeff3, coeff4, sphi, cphi, x52, d31, d41,
            d32, d42; 

    xm = x[0] - p->D * cos(p->omega * t);
    ym = x[1] - p->D * sin(p->omega * t);
    xm2 = xm * xm;
    ym2 = ym * ym;
    x12 = x[0] * x[0];
    x22 = x[1] * x[1];
    rE = sqrt(x12 + x22);
    rM = sqrt(xm2 + ym2);
    rE2 = rE * rE;
    rM2 = rM * rM;
    rE3 = rE2 * rE;
    rM3 = rM2 * rM;    
    coeff1 = p->gE / rE3;
    coeff2 = p->gM / rM3;
    coeff3 = coeff1 / rE2;
    coeff4 = coeff2 / rM2;
    x52 = x[4] * x[4];
    sphi = sin(u[1]);
    cphi = cos(u[1]);
    
	d31 = coeff3 * (rE2 - 3 * x12) + coeff4 * (rM2 - 3 * xm2);
	d41 = - 3 * (coeff3 * x[0] * x[1] + coeff4 * xm * ym);
	d32 = d41;
	d42 = coeff3 * (rE2 - 3 * x22) + coeff4 * (rM2 - 3 * ym2);
    
    out[0] = x[2];
    out[1] = x[3];
    out[2] = - coeff1 * x[0] - coeff2 * xm + u[0] * cos(u[1]) / x[4];
    out[3] = - coeff1 * x[1] - coeff2 * ym + u[0] * sin(u[1]) / x[4];
    out[4] = p->C1 * u[0];
    // rownania sprzezone
	out[5] = d31 * x[7] + d41 * x[8];
	out[6] = d32 * x[7] + d42 * x[8];
	out[7] = - x[5];
	out[8] = - x[6];
	out[9] = u[0] * (x[7] * cphi + x[8] * sphi) / x52;
    // dH/du
    out[10] = (x[7] * cphi + x[8] * sphi) / x[4] + x[9] * p->C1;
	out[11] = u[0] * (- x[7] * sphi + x[8] * cphi) / x[4];
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
    
    plhs[0] = mxCreateDoubleMatrix(1, 12, mxREAL);
 
    t = mxGetScalar(prhs[0]);
    x = mxGetPr(prhs[1]);
    u = mxGetPr(prhs[2]);
    out = mxGetPr(plhs[0]);
    
    rhsPsi(out, t, x, u, &p);
}
