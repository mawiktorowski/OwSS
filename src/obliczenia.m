function [ var ] = obliczenia(h0, tau)
% funkcja wyliczajaca potrzebne dane takie jak h, t, cn itd.

dtau = diff(tau);
n = ceil(dtau./h0);
h = dtau./n;
h2 = h/2;
h3 = h/3;
h6 = h/6;
cn = cumsum([1 n]);
iter = cn(end);
ldtau = length(dtau);
tf = tau(end);

var = struct('h', h, 'h2', h2, 'h3', h3, 'h6', h6, 'cn', cn, 'iter', iter, 'ldtau', ldtau, 'tf', tf);
end
