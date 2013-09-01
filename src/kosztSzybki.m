function [ Q, kara ] = kosztSzybki(zd, param, var, rho)
% szybkie wyliczanie wskaznika jakosci
        
        x = zeros(1,5);
        us = zeros(2,1);
        u = zd(3:end);
        
        x(1) = param.rE * cos(zd(1)) - param.mu;
        x(2) = param.rE * sin(zd(1));
        x(3) = - (param.VE + zd(2)) * sin(zd(1));
        x(4) = (param.VE + zd(2)) * cos(zd(1));
        x(5) = param.m0 * exp(param.C2 * zd(2));
        
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