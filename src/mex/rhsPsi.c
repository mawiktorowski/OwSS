#include "mex.h" /* Always include this */
#include <math.h>

struct param{
    double gE, gM, D3, C1;
};

void rhsPsi(double out[], double x[], double u[], struct param *p)
{
    double xm, ym, xm2, ym2, x32, x42, rE, rM, rE2, rM2, rE3, rM3,
            coeff1, coeff2, coeff3, coeff4, coeff5, sphi, cphi, x92,
            d71, d81, d72, d82, d73, d83, d74, d84; 

    xm = x[2] - x[0];
    ym = x[3] - x[1];
    xm2 = xm * xm;
    ym2 = ym * ym;
    x32 = x[2] * x[2];
    x42 = x[3] * x[3];
    rE = sqrt(x32 + x42);
    rM = sqrt(xm2 + ym2);
    rE2 = rE * rE;
    rM2 = rM * rM;
    rE3 = rE2 * rE;
    rM3 = rM2 * rM;
    coeff1 = p->gE / p->D3;    
    coeff2 = p->gE / rE3;
    coeff3 = p->gM / rM3;
    coeff4 = coeff2 / rE2;
    coeff5 = coeff3 / rM2;
    x92 = x[8] * x[8];
    sphi = sin(u[1]);
    cphi = cos(u[1]);
    
	d71 = - coeff5 * (rM2 - 3 * xm2);
	d81 = coeff5 * 3 * xm * ym;
	d72 = d81;
	d82 = - coeff5 * (rM2 - 3 * ym2);
	d73 = coeff4 * (rE2 - 3 * x32) + coeff5 * (rM2 - 3 * xm2);
	d83 = - 3 * (coeff4 * x[2] * x[3] + coeff5 * xm * ym);
	d74 = d83;
	d84 = coeff4 * (rE2 - 3 * x42) + coeff5 * (rM2 - 3 * ym2);
    

        out[0] = x[4];
    out[1] = x[5];
    out[2] = x[6];
    out[3] = x[7];
    out[4] = - coeff1 * x[0];
    out[5] = - coeff1 * x[1];
    out[6] = - coeff2 * x[2] - coeff3 * xm + u[0] * cphi / x[8];
    out[7] = - coeff2 * x[3] - coeff3 * ym + u[0] * sphi / x[8];
    out[8] = p->C1 * u[0];
    // rownania sprzezone
    out[9] = coeff1 * x[13] + d71 * x[15] + d81 * x[16];
	out[10] = coeff1 * x[14] + d72 * x[15] + d82 * x[16];
	out[11] = d73 * x[15] + d83 * x[16];
	out[12] = d74 * x[15] + d84 * x[16];
	out[13] = - x[9];
	out[14] = - x[10];
	out[15] = - x[11];
	out[16] = - x[12];
	out[17] = u[0] * (x[15] * cphi + x[16] * sphi) / x92;
    // dH/du
    out[18] = (x[15] * cphi + x[16] * sphi) / x[8] + x[17] * p->C1;
	out[19] = u[0] * (- x[15] * sphi + x[16] * cphi) / x[8];
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
    
    plhs[0] = mxCreateDoubleMatrix(1, 20, mxREAL);
    
    x = mxGetPr(prhs[0]);
    u = mxGetPr(prhs[1]);
    out = mxGetPr(plhs[0]);
    
    rhsPsi(out, x, u, &p);
}
