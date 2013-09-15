function [ zd, Q, kara ] = BFGS(zd, tau, param, rho, opts)
% algorytm BFGS z wbudowanym poszukiwaniem na kierunku

% --- Krok 1 ------------------------------------
var = obliczenia(param.h0, tau);
lzd = length(zd);
lu = (lzd - 2)/2;
urend = lu + 2;
phistart = urend + 1;
% epK = opts.epK0; % matlab twierdzi ze linia nigdzie nie uzywana

decmin = zeros(lzd,1);
decmin(1:2)          = param.ogr(1:2,1);
decmin(3:urend)      =  param.ogr(3,1);
decmin(phistart:lzd) = param.ogr(4,1);

decmax = zeros(lzd,1);
decmax(1:2)          = param.ogr(1:2,2);
decmax(3:urend)      =  param.ogr(3,2);
decmax(phistart:lzd) = param.ogr(4,2);

R = true;

[ Q, ~ ] = kosztSzybki(zd, param, var, rho);

for iter = 1:opts.MAX_ITER
    % --- Krok 2 ------------------------------------
    grad = solverSzybki(zd, param, var, rho);
    grad( (zd < decmin + opts.epOgr) & (grad > 0)) = 0;
    grad( (zd > decmax - opts.epOgr) & (grad < 0)) = 0;
    norma = norm(grad);
    
    if norma <= opts.ep0
        %fprintf(['KONIEC - MALA NORMA GRADIENTU\nITERACJA: ', num2str(iter) '\n']);
        break
    end
    
    % --- Krok 3 ------------------------------------
    if R
        W = eye(length(zd));
        epK = opts.epK0;
    else
        r = grad - gradOld;
        s = zd - zdOld;
        W = W + ( r * r' )/( s' * r ) - ( W * s * s' * W)/(s' * W * s);
        epK = opts.epK1;
    end
    
    % --- Krok 4 ------------------------------------
    d = mldivide(W,-grad);
    
    % --- Krok 5 ------------------------------------
    if d'*grad > -max(opts.ep1, opts.ep2 * norma^2)
        R = true;
        continue;
    end
    
    % --- Krok 6 ------------------------------------
    zdOld = zd;
    gradOld = grad;
    % --- Line Search -------------------------------
    lambda = 1;
    while(lambda > epK)
        zdNew = zd + lambda * d;
        zdNew = max(zdNew, decmin);
        zdNew = min(zdNew, decmax);
        
        [ QNew, kara ] = kosztSzybki(zdNew, param, var, rho);
        if(QNew < Q)
            zd = zdNew;
            Q = QNew;
            break;
        else
            lambda = lambda / 2; % kontrakcja
        end
    end
    % --- Krok 7 ------------------------------------
    if max(abs(zdOld - zd))
        R = false;
        %disp([Q norma]);
        continue;
    else
        if R
            %disp(iter);
            %fprintf('KONIEC - KIERUNEK NAJSZYBSZEGO SPADKU I BRAK POPRAWY\n');
            break;
        else
            %fprintf('ODNOWA ALGORYTMU\n');
            R = true;
        end
    end
end
end
