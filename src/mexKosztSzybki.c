#include "mex.h"
#include <math.h>

struct param{
    double mu, rmu, rE, VE, rM2, VM2, m0, mr, C1, C2, K1, K2;
};

struct var{
    double ldtau, *h, *h2, *h3, *h6, *cn, tf;
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

// function [ Q, kara ] = kosztSzybki(zd, param, var, rho)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    double *ptrZd, x[5], rho, u[2], hj, h2j, h3j, h6j, dx1[5], dx2[5],
            dx3[5], dx4[5], farg[5], xw, xw2, yw2, uw2, vw2, rsk1, rsk2,
            rsk3, k1, k2, k3, k4, *kara, *Q;
    int i, j, k, jj;
    struct param p;
    struct var v;
    
    ptrZd = mxGetPr(prhs[0]);
    
    p.mu = mxGetScalar(mxGetField(prhs[1], 0, "mu"));
    p.rmu = mxGetScalar(mxGetField(prhs[1], 0, "restmu"));
    p.rE = mxGetScalar(mxGetField(prhs[1], 0, "rE"));
    p.VE = mxGetScalar(mxGetField(prhs[1], 0, "VE"));
    p.rM2 = mxGetScalar(mxGetField(prhs[1], 0, "rM2"));
    p.VM2 = mxGetScalar(mxGetField(prhs[1], 0, "VM2"));
    p.m0 = mxGetScalar(mxGetField(prhs[1], 0, "m0"));
    p.mr = mxGetScalar(mxGetField(prhs[1], 0, "mr"));
    p.C1 = mxGetScalar(mxGetField(prhs[1], 0, "C1"));
    p.C2 = mxGetScalar(mxGetField(prhs[1], 0, "C2"));
    p.K1 = mxGetScalar(mxGetField(prhs[1], 0, "K1"));
    p.K2 = mxGetScalar(mxGetField(prhs[1], 0, "K2"));
    
    // struktura v //
    
    v.ldtau = mxGetScalar(mxGetField(prhs[2], 0, "ldtau"));
    v.h = mxGetPr(mxGetField(prhs[2], 0, "h"));
    v.h2 = mxGetPr(mxGetField(prhs[2], 0, "h2"));
    v.h3 = mxGetPr(mxGetField(prhs[2], 0, "h3"));
    v.h6 = mxGetPr(mxGetField(prhs[2], 0, "h6"));
    v.cn = mxGetPr(mxGetField(prhs[2], 0, "cn"));
    v.tf = mxGetScalar(mxGetField(prhs[2], 0, "tf"));
    
    rho = mxGetScalar(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    Q = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    kara = mxGetPr(plhs[1]);
    
    x[0] = p.rE * cos(ptrZd[0]) - p.mu;
    x[1] = p.rE * sin(ptrZd[0]);
    x[2] = - (p.VE + ptrZd[1]) * sin(ptrZd[0]);
    x[3] = (p.VE + ptrZd[1]) * cos(ptrZd[0]);
    x[4] = p.m0 * exp(p.C2 * ptrZd[1]);
    
    for(j=0; j < v.ldtau; j++){
        jj = v.ldtau + j;
        u[0] = ptrZd[j+2];
        u[1] = ptrZd[jj+2];
        hj =  v.h[j];
        h2j = v.h2[j];
        h3j = v.h3[j];
        h6j = v.h6[j];
        //printf("%f\n", v.h[j]);
        for(i = v.cn[j]; i < v.cn[j+1]; i++){
            // calosc petli do przerobienia potem
            rhs(dx1, x, u, &p);
            for(k=0; k < 5; k++) farg[k] = x[k] + h2j * dx1[k];
            rhs(dx2, farg, u, &p);
            for(k=0; k < 5; k++) farg[k] = x[k] + h2j * dx2[k];
            rhs(dx3, farg, u, &p);
            for(k=0; k < 5; k++) farg[k] = x[k] + hj * dx3[k];
            rhs(dx4, farg, u, &p);
            for(k=0; k < 5; k++) x[k] = x[k] + h3j * (dx2[k] + dx3[k]) + h6j * (dx1[k] + dx4[k]);
        }
    }
    
    for(k=0; k < 5; k++) printf("%f ", x[k]);
    
    xw = x[0] - p.rmu;
    
    xw2 = xw * xw;
    yw2 = x[1] * x[1];
    uw2 = x[2] * x[2];
    vw2 = x[3] * x[3];
    
    rsk1 = xw2 + yw2 - p.rM2;
    rsk2 = uw2 + vw2 - p.VM2;
    rsk3 = xw * x[2] + x[1] * x[3];
    
    k1 = 0.25 * rsk1 * rsk1;
    k2 = 0.25 * rsk2 * rsk2;
    k3 = 0.5 * rsk3 * rsk3;
    
    if (x[4] < p.mr) k4 = 0.5 * pow((x[4] - p.mr), 2);
    else k4 = 0;
    
    *kara = k1 + k2 + k3 + k4;
    *Q = - p.K1 * x[4] + p.K2 * v.tf + rho * *kara;
    
}