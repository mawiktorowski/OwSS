function [ dQ, H ] = testPsi(x0, u, param, var, ep, rho)
% sprawdzenie poprawnosci rozwiazywania rownan sprzezonych

dQ = zeros(5,1);

Q0 = testPsiKoszt(x0);
for m=1:5
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
        t = 0;
        
        for j = 1:var.ldtau
            jj = var.ldtau + j;
            us(1) = u(j);
            us(2) = u(jj);
            hj = var.h(j);
            h2j = var.h2(j);
            h3j = var.h3(j);
            h6j = var.h6(j);
            for i = var.cn(j):var.cn(j+1)-1
                dx1 = rhs(t, x, us, param);
                targ = t + h2j;
                farg = x + h2j * dx1;
                dx2 = rhs(targ, farg, us, param);
                farg = x + h2j * dx2;
                dx3 = rhs(targ, farg, us, param);
                targ = t + hj;
                farg = x + hj * dx3;
                dx4 = rhs(targ, farg, us, param);
                x = x + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
                t = t + hj;
            end
        end

        xm = x(1) -               param.D * cos(param.omega * t);
        ym = x(2) -               param.D * sin(param.omega * t);
        um = x(3) + param.omega * param.D * sin(param.omega * t);
        vm = x(4) - param.omega * param.D * cos(param.omega * t);
        
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
        if x(5) < param.mr
            k4 = 0.5 * (x(5) - param.mr)^2;
        else
            k4 = 0;
        end
        
        kara = k1 + k2 + k3 + k4;
        Q = - param.K1 * x(5) + param.K2 * var.tf + rho * kara;
        
    end

    function H = testPsiSolver(arg)
        % solver Psi dla sprawdzenia
        
        us = zeros(2,1);
        x = arg;
        t = 0;
        
        for j = 1:var.ldtau
            jj = var.ldtau + j;
            us(1) = u(j);
            us(2) = u(jj);
            hj = var.h(j);
            h2j = var.h2(j);
            h3j = var.h3(j);
            h6j = var.h6(j);
            for i = var.cn(j):var.cn(j+1)-1
                dx1 = rhs(t, x, us, param);
                targ = t + h2j;
                farg = x + h2j * dx1;
                dx2 = rhs(targ, farg, us, param);
                farg = x + h2j * dx2;
                dx3 = rhs(targ, farg, us, param);
                farg = x + hj * dx3;
                targ = t + hj;
                dx4 = rhs(targ, farg, us, param);
                x = x + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
                t = t + hj;
            end
        end
        
        xm = x(1) -               param.D * cos(param.omega * t);
        ym = x(2) -               param.D * sin(param.omega * t);
        um = x(3) + param.omega * param.D * sin(param.omega * t);
        vm = x(4) - param.omega * param.D * cos(param.omega * t);
        
        xm2 = xm * xm;
        ym2 = ym * ym;
        um2 = um * um;
        vm2 = vm * vm;
        
        beta1 = xm2 + ym2 - param.rM2;
        beta2 = um2 + vm2 - param.VM2;
        beta3 = xm * um + ym * vm;
        
        if x(5) <= param.mr
            diffK4 = x(5) - param.mr;
        else
            diffK4 = 0;
        end
        
        psi = zeros(1,5);
        
        psi(1) =  - rho * ( beta1 * xm + beta3 * um);
        psi(2) =  - rho * ( beta1 * ym + beta3 * vm);
        psi(3) =  - rho * ( beta2 * um + beta3 * xm);
        psi(4) =  - rho * ( beta2 * vm + beta3 * ym);
        psi(5) = param.K1 - rho * diffK4;
        
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
                dx1 = rhsPsi(t, x, us, param);
                targ = t - h2j;
                farg = x - h2j * dx1;
                dx2 = rhsPsi(targ, farg, us, param);
                farg = x - h2j * dx2;
                dx3 = rhsPsi(targ, farg, us, param);
                targ = t - hj;
                farg = x - hj * dx3;
                dx4 = rhsPsi(targ, farg, us, param);
                x = x - h3j * (dx2 + dx3) - h6j * (dx1 + dx4);
                t = t - hj;
            end
        end
        
        H = -x(6:10)';
    end

end
