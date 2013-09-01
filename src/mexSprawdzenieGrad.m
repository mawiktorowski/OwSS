function [ dQ, grad ] = sprawdzenieGrad(zd, param, var, ep, rho)

dQ = zeros(length(zd), 1);

Q0 = sprawdzenieGradKoszt(zd);

for m=1:length(zd)
    zdtmp = zd;
    zdtmp(m) = zd(m)+ep;
    dQ(m) = sprawdzenieGradKoszt(zdtmp);
end

dQ = (dQ - Q0) ./ ep;

grad = sprawdzenieGradSolver(zd);

    function Q = sprawdzenieGradKoszt(arg)
        % szybkie wyliczanie kosztu dla sprawdzenia
        
        x = zeros(1,5);
        us = zeros(2,1);
        u = arg(3:end);
        
        x(1) = param.rE * cos(arg(1)) - param.mu;
        x(2) = param.rE * sin(arg(1));
        x(3) = - (param.VE + arg(2)) * sin(arg(1));
        x(4) = (param.VE + arg(2)) * cos(arg(1));
        x(5) = param.m0 * exp(param.C2 * arg(2));
        
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
                farg = dgemm(x, h2j, dx1);
                % farg = x + h2j * dx1;
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

    function grad = sprawdzenieGradSolver(arg)
        % solver gradientu dla sprawdzenia
        
        x = zeros(1,5);
        gradU = zeros(2 * var.ldtau,1);
        us = zeros(2,1);
        u = arg(3:end);
        
        x(1) = param.rE * cos(arg(1)) - param.mu;
        x(2) = param.rE * sin(arg(1));
        % zd(3) potem do zmiany na zd(2)
        x(3) = - (param.VE + arg(2)) * sin(arg(1));
        x(4) = (param.VE + arg(2)) * cos(arg(1));
        x(5) = param.m0 * exp(param.C2 * arg(2));
        
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
            gradU(j) = x(11);
            gradU(jj) = x(12);
            x(11:12) = [0 0]; % gamma(t_i+1) = 0
        end
        
        dx1theta = - param.rE * sin(arg(1));
        dx2theta = param.rE * cos(arg(1));
        dx3theta = - (arg(2) + param.VE) * cos(arg(1));
        dx4theta = - (arg(2) + param.VE) * sin(arg(1));
        dx3dV = - sin(arg(1));
        dx4dV = cos(arg(1));
        dx5dV = param.m0 * param.C2 * exp(param.C2 * arg(2));
        
        grad = zeros(length(arg),1);
        grad(1) = - x(6) * dx1theta - x(7) * dx2theta - x(8) * dx3theta - x(9) * dx4theta;
        grad(2) = - x(8) * dx3dV - x(9) * dx4dV - x(10) * dx5dV;
        grad(3:end) = gradU;
        
    end

end