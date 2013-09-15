function [ dQ, H ] = testPsi(x0, u, param, var, ep, rho)
% sprawdzenie poprawnosci rozwiazywania rownan sprzezonych

dQ = zeros(9,1);

Q0 = testPsiKoszt(x0);
for m=1:9
    x0tmp = x0;
    x0tmp(m) = x0(m)+ep;
    dQ(m) = testPsiKoszt(x0tmp);
end

dQ = (dQ - Q0) ./ ep;

H = testPsiSolver(x0);

    function Q = testPsiKoszt(arg)
        % szybkie wyliczanie kosztu dla sprawdzenia
        
        us = zeros(2,1);
        x = arg;
        
        for j = 1:var.ldtau
            jj = var.ldtau + j;
            us(1) = u(j);
            us(2) = u(jj);
            hj = var.h(j);
            h2j = var.h2(j);
            h3j = var.h3(j);
            h6j = var.h6(j);
            for i = var.cn(j):var.cn(j+1)-1
                dx1 = rhs(x, us, param);
                farg = x + h2j * dx1;
                dx2 = rhs(farg, us, param);
                farg = x + h2j * dx2;
                dx3 = rhs(farg, us, param);
                farg = x + hj * dx3;
                dx4 = rhs(farg, us, param);
                x = x + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
            end
        end

        xm = x(3) - x(1);
        ym = x(4) - x(2);
        um = x(7) - x(5);
        vm = x(8) - x(6);
        
        xm2 = xm * xm;
        ym2 = ym * ym;
        um2 = um * um;
        vm2 = vm * vm;
        
        beta1 = xm2 + ym2 - param.rM2;
        beta2 = um2 + vm2 - param.VM2;
        beta3 = xm * um + ym * vm;
        
        k1 = 0.25 * beta1 * beta1;
        k2 = 0.25 * beta2 * beta2;
        k3 = 0.5 * beta3 * beta3;
        if x(9) < param.mr
            k4 = 0.5 * (x(9) - param.mr)^2;
        else
            k4 = 0;
        end
        
        kara = k1 + k2 + k3 + k4;
        Q = - param.K1 * x(9) + param.K2 * var.tf + rho * kara;
        
    end

    function H = testPsiSolver(arg)
        % solver Psi dla sprawdzenia
        
        us = zeros(2,1);
        x = arg;
        
        for j = 1:var.ldtau
            jj = var.ldtau + j;
            us(1) = u(j);
            us(2) = u(jj);
            hj = var.h(j);
            h2j = var.h2(j);
            h3j = var.h3(j);
            h6j = var.h6(j);
            for i = var.cn(j):var.cn(j+1)-1
                dx1 = rhs(x, us, param);
                farg = x + h2j * dx1;
                dx2 = rhs(farg, us, param);
                farg = x + h2j * dx2;
                dx3 = rhs(farg, us, param);
                farg = x + hj * dx3;
                dx4 = rhs(farg, us, param);
                x = x + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
            end
        end
        
        xm = x(3) - x(1);
        ym = x(4) - x(2);
        um = x(7) - x(5);
        vm = x(8) - x(6);
        
        xm2 = xm * xm;
        ym2 = ym * ym;
        um2 = um * um;
        vm2 = vm * vm;
        
        beta1 = xm2 + ym2 - param.rM2;
        beta2 = um2 + vm2 - param.VM2;
        beta3 = xm * um + ym * vm;
        
        if x(9) <= param.mr
            diffK4 = x(9) - param.mr;
        else
            diffK4 = 0;
        end
        
        psi = zeros(1,9);
        
        psi(1) =   rho * ( beta1 * xm + beta3 * um);
        psi(2) =   rho * ( beta1 * ym + beta3 * vm);
        psi(3) = - rho * ( beta1 * xm + beta3 * um);
        psi(4) = - rho * ( beta1 * ym + beta3 * vm);
        psi(5) =   rho * ( beta2 * um + beta3 * xm);
        psi(6) =   rho * ( beta2 * vm + beta3 * ym);
        psi(7) = - rho * ( beta2 * um + beta3 * xm);
        psi(8) = - rho * ( beta2 * vm + beta3 * ym);
        psi(9) = param.K1 - rho * diffK4;
        
        x = [x psi 0 0];
        
        for j = var.ldtau:-1:1
            jj = var.ldtau + j;
            us(1) = u(j);
            us(2) = u(jj);
            hj = var.h(j);
            h2j = var.h2(j);
            h3j = var.h3(j);
            h6j = var.h6(j);
            for i = (var.cn(j+1)-1):-1:var.cn(j)
                dx1 = rhsPsi(x, us, param);
                farg = x - h2j * dx1;
                dx2 = rhsPsi(farg, us, param);
                farg = x - h2j * dx2;
                dx3 = rhsPsi(farg, us, param);
                farg = x - hj * dx3;
                dx4 = rhsPsi(farg, us, param);
                x = x - h3j * (dx2 + dx3) - h6j * (dx1 + dx4);
            end
        end
        
        H = -x(10:18)';
    end

end
