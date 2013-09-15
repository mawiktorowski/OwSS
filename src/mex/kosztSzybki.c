#include "mex.h"
#include <math.h>

struct param{
    double gE, gM, D, D3, omega, rE, VE, rM2, VM2, m0, mr, C1, C2, K1, K2;
};

struct var{
    double ldtau, *h, *h2, *h3, *h6, *cn, tf;
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

// function [ Q, kara ] = kosztSzybki(zd, param, var, rho)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    double *ptrZd, x[9], rho, u[2], hj, h2j, h3j, h6j, dx1[9], dx2[9],
            dx3[9], dx4[9], farg[9], xm, ym, um, vm ,xm2, ym2, um2, vm2,
            beta1, beta2, beta3, k1, k2, k3, k4, *kara, *Q;
    int i, j, k, jj;
    struct param p;
    struct var v;
    
    ptrZd = mxGetPr(prhs[0]);
    
    p.gE = mxGetScalar(mxGetField(prhs[1], 0, "gE"));
    p.gM = mxGetScalar(mxGetField(prhs[1], 0, "gM"));
    p.D = mxGetScalar(mxGetField(prhs[1], 0, "D"));
    p.D3 = mxGetScalar(mxGetField(prhs[1], 0, "D3"));
    p.omega = mxGetScalar(mxGetField(prhs[1], 0, "omega"));
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
    
    x[0] = p.D;
    x[1] = 0;
    x[2] = p.rE * cos(ptrZd[0]);
    x[3] = p.rE * sin(ptrZd[0]);
    x[4] = 0;
    x[5] = p.omega * p.D;
    x[6] = - (p.VE + ptrZd[1]) * sin(ptrZd[0]);
    x[7] = (p.VE + ptrZd[1]) * cos(ptrZd[0]);
    x[8] = p.m0 * exp(p.C2 * ptrZd[1]);
    
    for(j=0; j < v.ldtau; j++){
        jj = v.ldtau + j;
        u[0] = ptrZd[j+2];
        u[1] = ptrZd[jj+2];
        hj =  v.h[j];
        h2j = v.h2[j];
        h3j = v.h3[j];
        h6j = v.h6[j];
        for(i = v.cn[j]; i < v.cn[j+1]; i++){
            // calosc petli do przerobienia potem
            rhs(dx1, x, u, &p);
            for(k=0; k < 9; k++) farg[k] = x[k] + h2j * dx1[k];
            rhs(dx2, farg, u, &p);
            for(k=0; k < 9; k++) farg[k] = x[k] + h2j * dx2[k];
            rhs(dx3, farg, u, &p);
            for(k=0; k < 9; k++) farg[k] = x[k] + hj * dx3[k];
            rhs(dx4, farg, u, &p);
            for(k=0; k < 9; k++) x[k] = x[k] + h3j * (dx2[k] + dx3[k]) + h6j * (dx1[k] + dx4[k]);
        }
    }
    
    xm = x[2] - x[0];
    ym = x[3] - x[1];
    um = x[6] - x[4];
    vm = x[7] - x[5];

    xm2 = xm * xm;
    ym2 = ym * ym;
    um2 = um * um;
    vm2 = vm * vm;
    
    beta1 = xm2 + ym2 - p.rM2;
    beta2 = um2 + vm2 - p.VM2;
    beta3 = xm * um + ym * vm;
    
    k1 = 0.25 * beta1 * beta1;
    k2 = 0.25 * beta2 * beta2;
    k3 = 0.5 * beta3 * beta3;
    
    if (x[8] < p.mr) k4 = 0.5 * pow((x[8] - p.mr), 2);
    else k4 = 0;
    
    *kara = k1 + k2 + k3 + k4;
    *Q = - p.K1 * x[8] + p.K2 * v.tf + rho * *kara;
    
}