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

// function [ Q, grad ] = solverSzybki(zd, param, var, rho)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    double *ptrZd, x[9], rho, u[2], hj, h2j, h3j, h6j, dx1[9], dx2[9],
            dx3[9], dx4[9], farg[9], xm, ym, um, vm, xm2, ym2, um2, vm2,
            beta1, beta2, beta3, *grad, diffK4,
            psi[20], dx3theta, dx4theta, dx7theta, dx8theta, dx7dV, dx8dV,
            dx9dV, pdx1[20], pdx2[20], pdx3[20], pdx4[20], pfarg[20];
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
    
    plhs[0] = mxCreateDoubleMatrix((2 + 2 * v.ldtau), 1, mxREAL);
    grad = mxGetPr(plhs[0]);
    
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
    
    if (x[8] <= p.mr) diffK4 = x[8] - p.mr;
    else diffK4 = 0;
    
    for(k=0; k < 9; k++) psi[k] = x[k];
    
    psi[9]  =   rho * ( beta1 * xm + beta3 * um);
    psi[10] =   rho * ( beta1 * ym + beta3 * vm);
    psi[11]  = - psi[9] ;
    psi[12]  = - psi[10];
    psi[13]  =   rho * ( beta2 * um + beta3 * xm);
    psi[14]  =   rho * ( beta2 * vm + beta3 * ym);
    psi[15]  = - psi[13];
    psi[16]  = - psi[14];
    psi[17]  = p.K1 - rho * diffK4;
    psi[18] = 0;
    psi[19] = 0;
    
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
            rhsPsi(pdx1, psi, u, &p);
            for(k=0; k < 20; k++) pfarg[k] = psi[k] - h2j * pdx1[k];
            rhsPsi(pdx2, pfarg, u, &p);
            for(k=0; k < 20; k++) pfarg[k] = psi[k] - h2j * pdx2[k];
            rhsPsi(pdx3, pfarg, u, &p);
            for(k=0; k < 20; k++) pfarg[k] = psi[k] - hj * pdx3[k];
            rhsPsi(pdx4, pfarg, u, &p);
            for(k=0; k < 20; k++) psi[k] = psi[k] - h3j * (pdx2[k] + pdx3[k]) - h6j * (pdx1[k] + pdx4[k]);
        }
        grad[j+2] = psi[18];
        psi[18] = 0;
        grad[jj+2] = psi[19];
        psi[19] = 0;
    }
    
    dx3theta = - p.rE * sin(ptrZd[0]);
    dx4theta = p.rE * cos(ptrZd[0]);
    dx7theta = - (ptrZd[1] + p.VE) * cos(ptrZd[0]);
    dx8theta = - (ptrZd[1] + p.VE) * sin(ptrZd[0]);
    dx7dV = - sin(ptrZd[0]);
    dx8dV = cos(ptrZd[0]);
    dx9dV = p.m0 * p.C2 * exp(p.C2 * ptrZd[1]);

    grad[0] = - psi[11] * dx3theta - psi[12] * dx4theta - psi[15] * dx7theta - psi[16] * dx8theta;
    grad[1] = - psi[15] * dx7dV - psi[16] * dx8dV - psi[17] * dx9dV;  
}