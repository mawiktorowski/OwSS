function [ t, x, psi, grad ] = solver(x0, h0, tau, u)
% solver modelu rozwiazujacy rownania w przod i wstecz z zapamietywaniem

global restmu K2 rM VM mr

dtau = diff(tau);
n = ceil(dtau./h0);
h = dtau./n;
h2 = h/2;
h3 = h/3;
h6 = h/6;
cn = cumsum([1 n]);
iter = cn(end);

x = zeros(iter,length(x0));
t = zeros(iter, 1);
grad = zeros(2*length(dtau),1);
us = zeros(2,1);

x(1,:) = x0;

for j = 1:length(dtau)
    us(1) = u(j);
    us(2) = u(length(dtau) + j);   
    f = @(x,u) rhs(x,u);
    for i = cn(j):cn(j+1)-1
        dx1 = f(x(i,:), us);
        dx2 = f(x(i,:) + h2(j) * dx1, us);
        dx3 = f(x(i,:) + h2(j) * dx2, us);
        dx4 = f(x(i,:) + h(j) * dx3, us);
        x(i+1,:) = x(i,:) + h3(j) * (dx2 + dx3) + h6(j) * (dx1 + dx4);
        t(i+1) = t(i) + h(j);
    end
end

T = t(iter);
xT = x(iter,:);

beta1 = rM^2 - (xT(1) - restmu)^2 - xT(2)^2;
beta2 = VM^2 - (xT(3) - xT(2))^2 - (xT(4) + xT(1))^2;
beta3 = - (xT(1) - restmu) * (xT(3) - xT(2)) - xT(2) * (xT(4) + xT(1));

    function y = diffK4
        if xT(5) <= mr
            y = xT(5) - mr;
        else
            y = 0;
        end
    end

psiT(1) = K2 * (beta1 * (xT(1) - restmu) + beta2 * (xT(4) + xT(1)) + beta3 * xT(3));
psiT(2) = K2 * (beta1 * xT(2)            + beta2 * (xT(2) - xT(3)) + beta3 * (restmu + xT(4)));
psiT(3) = K2 * (                          beta2 * (xT(3) - xT(2)) + beta3 * (xT(1) - restmu));
psiT(4) = K2 * (                          beta2 * (xT(4) + xT(1)) + beta3 * xT(2));
psiT(5) = 1 - K2 * diffK4;

xBack = zeros(iter, 2 * length(xT) + 2);
xBack(iter,:) = [xT psiT 0 0];

for j = length(dtau):-1:1
    us(1) = u(j);
    us(2) = u(length(dtau) + j);       
    g = @(x,u) rhsPsi(x,u);
    for i = (cn(j+1)-1):-1:cn(j)
        dx1 = g(xBack(i+1,:), us);
        dx2 = g(xBack(i+1,:) - h2(j) * dx1, us);
        dx3 = g(xBack(i+1,:) - h2(j) * dx2, us);
        dx4 = g(xBack(i+1,:) - h(j) * dx3, us);
        xBack(i,:) = xBack(i+1,:) - h3(j) * (dx2 + dx3) - h6(j) * (dx1 + dx4);
    end
    grad(j,1) = xBack(i,11);
    grad(length(dtau)+j,1) = xBack(i,12);
    xBack(i,11:12) = [0 0]; % gamma(t_i+1) = 0
end

psi = xBack(:,length(xT) + 1:2 * length(xT));
%psi = xBack;
end