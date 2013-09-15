function [ dH, psi ] = testH(t, x0, psi0, u, param, ep)
% sprawdzenie poprawnosci wyliczania rownan sprzezonych

dH = zeros(5,1);
H0 = sum(psi0 .* rhs(t, x0, u, param));

for m=1:5
    x0tmp = x0;
    x0tmp(m) = x0(m)+ep;
    dH(m) = sum(psi0 .* rhs(t, x0tmp, u, param));
end

dH = (dH - H0) ./ ep;
psi = rhsPsi(t, [x0 psi0], u, param);
psi = - psi(6:10)';

end
