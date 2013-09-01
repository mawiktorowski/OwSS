function [ dQ, H ] = sprawdzeniePsi(x0, u, param, var, ep, rho)
% sprawdzenie poprawnosci rozwiazywania rownan sprzezonych

dQ = zeros(5,1);

Q0 = sprawdzeniePsiKoszt(x0);
for m=1:5
    x0tmp = x0;
    x0tmp(m) = x0(m)+ep;
    dQ(m) = sprawdzeniePsiKoszt(x0tmp);
end

dQ = (dQ - Q0) ./ ep;

H = sprawdzeniePsiSolver(x0);

    function Q = sprawdzeniePsiKoszt(arg)
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
        
        xw = x(1) - param.restmu;
        
        xw2 = xw * xw;
        yw2 = x(2) * x(2);
        uw2 = x(3) * x(3);
        vw2 = x(4) * x(4);
        
        rsk1 = xw2 + yw2 - param.rM2;
        rsk2 = uw2 + vw2 - param.VM2;
        rsk3 = xw * x(3) + x(2) * x(4);
        
        k1 = 0.25 * rsk1 * rsk1;
        k2 = 0.25 * rsk2 * rsk2;
        k3 = 0.5 * rsk3 * rsk3;
        if x(5) < param.mr
            k4 = 0.5 * (x(5) - param.mr)^2;
        else
            k4 = 0;
        end
        
        kara = k1 + k2 + k3 + k4;
        Q = - param.K1 * x(5) + param.K2 * var.tf + rho * kara;
        
    end

function H = sprawdzeniePsiSolver(arg)
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

        xw = x(1) - param.restmu;
        
        xw2 = xw * xw;
        yw2 = x(2) * x(2);
        uw2 = x(3) * x(3);
        vw2 = x(4) * x(4);
        
        beta1 = xw2 + yw2 - param.rM2;
        beta2 = uw2 + vw2 - param.VM2;
        beta3 = xw * x(3) + x(2) * x(4);
        
        if x(5) <= param.mr
            diffK4 = x(5) - param.mr;
        else
            diffK4 = 0;
        end 
        
        psi = zeros(1,5);

        psi(1) = - rho * (beta1 * xw   + beta3 * x(3));
        psi(2) = - rho * (beta1 * x(2) + beta3 * x(4));
        psi(3) = - rho * (beta2 * x(3) + beta3 * xw);
        psi(4) = - rho * (beta2 * x(4) + beta3 * x(2));
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

        H = -x(6:10)';
end

end
