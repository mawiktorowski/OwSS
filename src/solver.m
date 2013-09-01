function [ x, psi, grad ] = solver(zd, param, var, rho)
% solver modelu rozwiazujacy rownania w przod i wstecz z zapamietywaniem

x = zeros(var.iter,5);
us = zeros(2,1);
u = zd(3:end);

x(1,1) = param.rE * cos(zd(1)) - param.mu;
x(1,2) = param.rE * sin(zd(1));
x(1,3) = - (param.VE + zd(2)) * sin(zd(1));
x(1,4) = (param.VE + zd(2)) * cos(zd(1));
x(1,5) = param.m0 * exp(param.C2 * zd(2));

for j = 1:var.ldtau
    jj = var.ldtau + j;
    us(1) = u(j);
    us(2) = u(jj);
    hj = var.h(j);
    h2j = var.h2(j);
    h3j = var.h3(j);
    h6j = var.h6(j);
    for i = var.cn(j):var.cn(j+1)-1
        xc = x(i,:);
        dx1 = rhs(xc, us, param);
        farg = xc + h2j * dx1;
        dx2 = rhs(farg, us, param);
        farg = xc + h2j * dx2;
        dx3 = rhs(farg, us, param);
        farg = xc + hj * dx3;
        dx4 = rhs(farg, us, param);
        x(i+1,:) = xc + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
    end
end

xT = x(var.iter,:);

xw = xT(1) - param.restmu;

xw2 = xw * xw;
yw2 = xT(2) * xT(2);
uw2 = xT(3) * xT(3);
vw2 = xT(4) * xT(4);

beta1 = xw2 + yw2 - param.rM2;
beta2 = uw2 + vw2 - param.VM2;
beta3 = xw * xT(3) + xT(2) * xT(4);

if xT(5) <= param.mr
    diffK4 = xT(5) - param.mr;
else
    diffK4 = 0;
end

psiT = zeros(1,5);

psiT(1) = - rho * (beta1 * xw   + beta3 * xT(3));
psiT(2) = - rho * (beta1 * xT(2) + beta3 * xT(4));
psiT(3) = - rho * (beta2 * xT(3) + beta3 * xw);
psiT(4) = - rho * (beta2 * xT(4) + beta3 * xT(2));
psiT(5) = param.K1 - rho * diffK4;

x = [x zeros(var.iter, 7)];

x(var.iter, 6:12) = [psiT 0 0];

for j = var.ldtau:-1:1
    jj = var.ldtau + j;
    us(1) = u(j);
    us(2) = u(jj);
    hj = var.h(j);
    h2j = var.h2(j);
    h3j = var.h3(j);
    h6j = var.h6(j);
    for i = (var.cn(j+1)-1):-1:var.cn(j)
        xc = x(i+1,:);
        dx1 = rhsPsi(xc, us, param);
        farg = xc - h2j * dx1;
        dx2 = rhsPsi(farg, us, param);
        farg = xc - h2j * dx2;
        dx3 = rhsPsi(farg, us, param);
        farg = xc - hj * dx3;
        dx4 = rhsPsi(farg, us, param);
        x(i,:) = xc - h3j * (dx2 + dx3) - h6j * (dx1 + dx4);
    end
    gradU(j) = x(i,11);
    gradU(jj) = x(i,12);
    x(i, 11:12) = [0 0]; % gamma(t_i+1) = 0
end

dx1theta = - param.rE * sin(zd(1));
dx2theta = param.rE * cos(zd(1));
dx3theta = - (zd(2) + param.VE) * cos(zd(1));
dx4theta = - (zd(2) + param.VE) * sin(zd(1));
dx3dV = - sin(zd(1));
dx4dV = cos(zd(1));
dx5dV = param.m0 * param.C2 * exp(param.C2 * zd(2));

grad = zeros(length(zd),1);
grad(1) = - x(1,6) * dx1theta - x(1,7) * dx2theta - x(1,8) * dx3theta - x(1,9) * dx4theta;
grad(2) = - x(1,8) * dx3dV - x(1,9) * dx4dV - x(1,10) * dx5dV;
grad(3:end) = gradU;

psi = x(:,6:12);
x = x(:,1:5);
end
