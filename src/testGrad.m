function [ dQ, grad ] = testGrad(zd, param, var, ep, rho)

dQ = zeros(length(zd), 1);

Q0 = kosztSzybki(zd, param, var, rho);

for m=1:length(zd)
    zdtmp = zd;
    zdtmp(m) = zd(m)+ep;
    dQ(m) = kosztSzybki(zdtmp, param, var, rho);
end

dQ = (dQ - Q0) ./ ep;

grad = solverSzybki(zd, param, var, rho);

end
