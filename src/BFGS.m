function [ zd, Q, kara ] = BFGS(zd, tau, param, rho, opts)
% algorytm BFGS z wbudowanym poszukiwaniem na kierunku

var = obliczenia(param.h0, tau);
lzd = length(zd);
lu = (lzd - 2)/2;
urend = lu + 2;
phistart = urend + 1;
% epK = opts.epK0; % matlab twierdzi ze linia nigdzie nie uzywana
R = true;
iter = 1;

while(iter <= opts.MAX_ITER)
    % --- Krok 1 ------------------------------------
    [Q, grad] = solverSzybki(zd, param, var, rho);
    %rzutowanie gradientu na ograniczenia
    grad(zd(1:2) < (param.ogr(1:2,1) + opts.epOgr) & grad(1:2) > 0) = 0;
    grad(zd(1:2) > (param.ogr(1:2,2) - opts.epOgr) & grad(1:2) < 0) = 0;
    grad(zd(3:urend) < (param.ogr(3,1) + opts.epOgr) & grad(3:urend) > 0) = 0;
    grad(zd(3:urend) > (param.ogr(3,2) - opts.epOgr) & grad(3:urend) < 0) = 0;
    grad(zd(phistart:end) < (param.ogr(4,1) + opts.epOgr) & grad(phistart:end) > 0) = 0;
    grad(zd(phistart:end) > (param.ogr(4,2) - opts.epOgr) & grad(phistart:end) < 0) = 0;
    norma = norm(grad);
    
    % --- Krok 2 ------------------------------------
    if norma <= opts.ep0
        fprintf(['KONIEC - MALA NORMA GRADIENTU\nITERACJA: ', num2str(iter) '\n']);
        break
    end
    
    step3 = true;
    while step3
        % --- Krok 3 ------------------------------------
        if R
            W = eye(length(zd));
            epK = opts.epK0;
        else
            r = grad - gradOld;
            s = zd - zdOld;
            Ws = W * s;
            W = W + (r * r')/(s' * r) - (Ws * s' * W)/(s' * Ws);
            epK = opts.epK1;
        end
        
        d = W \ (-grad);
        
        % --- Krok 4 ------------------------------------
        if d' * grad > - max(opts.ep1, opts.ep2 * norma^2)
            fprintf('ODNOWA ALGORYTMU\n');
            R = true;
            iter = iter + 1;
            continue
        end
        
        % --- Krok 5 ------------------------------------
        zdOld = zd;
        gradOld = grad;
        
        % --- Krok 6 ------------------------------------
        % --- Line Search -------------------------------
        % rzutowanie wektora d na ograniczenia
        %dOld = d;
        d(1:2) = max(d(1:2), param.ogr(1:2,1) - zd(1:2));
        d(1:2) = min(d(1:2), param.ogr(1:2,2) - zd(1:2));
        d(3:urend) = max(d(3:urend), param.ogr(3,1) - zd(3:urend));
        d(3:urend) = min(d(3:urend), param.ogr(3,2) - zd(3:urend));
        d(phistart:end) = max(d(phistart:end), param.ogr(4,1) - zd(phistart:end));
        d(phistart:end) = min(d(phistart:end), param.ogr(4,2) - zd(phistart:end));
        
        %disp([dOld d ogr(:,1) -zd ogr(:,2) - zd])
        lambda = 1;
        while(lambda > epK)
            zdNew = zd + lambda * d;
            [ QNew, kara ] = kosztSzybki(zdNew, param, var, rho);
            if(QNew < Q)
                zd = zdNew;
                Q = QNew;
                break
            else
                lambda = lambda / 2; % kontrakcja
            end
        end
        
        % --- Krok 7 ------------------------------------
        if max(abs(zdOld - zd))
            R = false;
            disp([Q norm(grad)]);
            break
        else
            if R
                fprintf('KONIEC - KIERUNEK NAJSZYBSZEGO SPADKU I BRAK POPRAWY\n');
                %max(abs(zdOld - zd));
                step3 = false;
                break
            else
                R = true;
                iter = iter + 1;
            end
        end
    end
    
    if ~step3
        break
    end
    iter = iter + 1;
    
end
end