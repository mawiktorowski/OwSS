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

// function [ Q, grad ] = solverSzybki(zd, param, var, rho)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    double *ptrZd, x[5], t, rho, u[2], hj, h2j, h3j, h6j, dx1[5], dx2[5],
            dx3[5], dx4[5], farg[5], targ, xm, ym, um, vm, xm2, ym2,
            um2, vm2, beta1, beta2, beta3, *grad, diffK4,
            psi[12], dx1theta, dx2theta, dx3theta, dx4theta, dx3dV, dx4dV,
            dx5dV, pdx1[12], pdx2[12], pdx3[12], pdx4[12], pfarg[12];
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
    
    plhs[0] = mxCreateDoubleMatrix((2 + 2 * v.ldtau), 1, mxREAL);
    grad = mxGetPr(plhs[0]);
    
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
    
    if (x[4] <= p.mr) diffK4 = x[4] - p.mr;
    else diffK4 = 0;
    
    for(k=0; k < 5; k++) psi[k] = x[k];
    
    psi[5]  = - rho * ( beta1 * xm + beta3 * um);
    psi[6]  = - rho * ( beta1 * ym + beta3 * vm);
    psi[7]  = - rho * ( beta2 * um + beta3 * xm);
    psi[8]  = - rho * ( beta2 * vm + beta3 * ym);
    psi[9]  = p.K1 - rho * diffK4;
    psi[10] = 0;
    psi[11] = 0;
    
      for(j=v.ldtau - 1; j >= 0; j--){
        jj = v.ldtau + j;
        u[0] = ptrZd[j+2];
        u[1] = ptrZd[jj+2];
        hj =  v.h[j];
        h2j = v.h2[j];
        h3j = v.h3[j];
        h6j = v.h6[j];
        for(i = v.cn[j+1] - 1; i >= v.cn[j]; i--){
            // calosc petli do przerobienia potem
            rhsPsi(pdx1, t, psi, u, &p);
            for(k=0; k < 12; k++) pfarg[k] = psi[k] - h2j * pdx1[k];
            targ = t - h2j;
            rhsPsi(pdx2, targ, pfarg, u, &p);
            for(k=0; k < 12; k++) pfarg[k] = psi[k] - h2j * pdx2[k];
            rhsPsi(pdx3, targ, pfarg, u, &p);
            for(k=0; k < 12; k++) pfarg[k] = psi[k] - hj * pdx3[k];
            targ = t - hj;
            rhsPsi(pdx4, targ, pfarg, u, &p);
            for(k=0; k < 12; k++) psi[k] = psi[k] - h3j * (pdx2[k] + pdx3[k]) - h6j * (pdx1[k] + pdx4[k]);
            t = t - hj;
        }
        grad[j+2] = psi[10];
        psi[10] = 0;
        grad[jj+2] = psi[11];
        psi[11] = 0;
    }
    
    dx1theta = - p.rE * sin(ptrZd[0]);
    dx2theta = p.rE * cos(ptrZd[0]);
    dx3theta = - (ptrZd[1] + p.VE) * cos(ptrZd[0]);
    dx4theta = - (ptrZd[1] + p.VE) * sin(ptrZd[0]);
    dx3dV = - sin(ptrZd[0]);
    dx4dV = cos(ptrZd[0]);
    dx5dV = p.m0 * p.C2 * exp(p.C2 * ptrZd[1]);

    grad[0] = - psi[5] * dx1theta - psi[6] * dx2theta - psi[7] * dx3theta - psi[8] * dx4theta;
    grad[1] = - psi[7] * dx3dV - psi[8] * dx4dV - psi[9] * dx5dV;  
}