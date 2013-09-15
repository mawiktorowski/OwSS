function [ Q, zd, tau, dop ] = funkcjaKary(zd, tau, param, opts)
% funkcja kary

dop = false;
rho = 1;

while rho < opts.k
    [ zdopt, ~, kara ] = BFGS(zd, tau, param, rho, opts);
    if (kara < opts.epKara)
        break
    end
    zd = zdopt;
    rho = rho * opts.c;
end

var = obliczenia(param.h0, tau);

[ Q, ~ ] = kosztSzybki(zd, param, var, 0);

if (kara < opts.epKara)
    dop = true;
    fprintf(['CEL OSI¥GNIÊTY     Q=', num2str(Q), ' T=', num2str(tau(end)) '\n']);
else
    fprintf(['CEL NIE OSI¥GNIÊTY Q=', num2str(Q), ' T=', num2str(tau(end)) '\n']);
end
end
