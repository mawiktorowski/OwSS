#include "mex.h"
#include <math.h>

struct param{
    double gE, gM, D, omega, rE, VE, rM2, VM2, m0, mr, C1, C2, K1, K2;
};

struct var{
    double ldtau, *h, *h2, *h3, *h6, *cn, tf;
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

// function [ Q, kara ] = kosztSzybki(zd, param, var, rho)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    double *ptrZd, x[5], t, rho, u[2], hj, h2j, h3j, h6j, dx1[5], dx2[5],
            dx3[5], dx4[5], farg[5], targ, xm, ym, um, vm ,xm2, ym2, um2, vm2,
            beta1, beta2, beta3, k1, k2, k3, k4, *kara, *Q;
    int i, j, k, jj;
    struct param p;
    struct var v;
    
    ptrZd = mxGetPr(prhs[0]);
    
    p.gE = mxGetScalar(mxGetField(prhs[1], 0, "gE"));
    p.gM = mxGetScalar(mxGetField(prhs[1], 0, "gM"));
    p.D = mxGetScalar(mxGetField(prhs[1], 0, "D"));
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
    
    x[0] = p.rE * cos(ptrZd[0]);
    x[1] = p.rE * sin(ptrZd[0]);
    x[2] = - (p.VE + ptrZd[1]) * sin(ptrZd[0]);
    x[3] = (p.VE + ptrZd[1]) * cos(ptrZd[0]);
    x[4] = p.m0 * exp(p.C2 * ptrZd[1]);

    t = 0;
    
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
            rhs(dx1, t, x, u, &p);
            for(k=0; k < 5; k++) farg[k] = x[k] + h2j * dx1[k];
            targ = t + h2j;
            rhs(dx2, targ, farg, u, &p);
            for(k=0; k < 5; k++) farg[k] = x[k] + h2j * dx2[k];
            rhs(dx3, targ, farg, u, &p);
            for(k=0; k < 5; k++) farg[k] = x[k] + hj * dx3[k];
            targ = t + hj;
            rhs(dx4, targ, farg, u, &p);
            for(k=0; k < 5; k++) x[k] = x[k] + h3j * (dx2[k] + dx3[k]) + h6j * (dx1[k] + dx4[k]);
            t = t + hj;
        }
    }
    
    xm = x[0] -           p.D * cos(p.omega * t);
    ym = x[1] -           p.D * sin(p.omega * t);
    um = x[2] + p.omega * p.D * sin(p.omega * t);
    vm = x[3] - p.omega * p.D * cos(p.omega * t);
    
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
    
    if (x[4] < p.mr) k4 = 0.5 * pow((x[4] - p.mr), 2);
    else k4 = 0;
    
    *kara = k1 + k2 + k3 + k4;
    *Q = - p.K1 * x[4] + p.K2 * v.tf + rho * *kara;
    
}