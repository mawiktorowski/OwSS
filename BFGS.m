function [ zd, Q ] = BFGS(x0, h0, tau, zd, ogr)
% algorytm BFGS z wbudowanym poszukiwaniem na kierunku

global ep0 ep1 ep2 epK0 epK1 MAX_ITER

epK = epK0;
R = true;
iter = 1;

f = @(x) solverSzybki(x0, h0, tau, x); % wyliczanie kosztu
g = @(x) solverSzybki(x0, h0, tau, x); % wyliczanie kosztu i gradientu

while(iter <= MAX_ITER)
    % --- Krok 1 ------------------------------------
    [Q, grad] = g(zd);
    %rzutowanie gradientu na ograniczenia
    grad(zd == ogr(:,1) & grad > 0) = 0;
    grad(zd == ogr(:,2) & grad < 0) = 0;
    
    % --- Krok 2 ------------------------------------
    if norm(grad) <= ep0
        fprintf(['KONIEC - MALA NORMA GRADIENTU\nITERACJA: ', num2str(iter) '\n']);
        break
    end

    step3 = true;
    while step3
        % --- Krok 3 ------------------------------------
        if R
            W = eye(length(zd));
            epK = epK0;
        else
            r = grad - gradOld;
            s = zd - zdOld;
            Ws = W * s;
            W = W + (r * r')/(s' * r) - (Ws * s' * W)/(s' * Ws);
            epK = epK1;
        end
        
        d = W \ (-grad);
        
        % --- Krok 4 ------------------------------------
        if d' * grad > - max(ep1, ep2 * norm(grad)^2)
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
        d = max(d, ogr(:,1) - zd);
        d = min(d, ogr(:,2) - zd);
        lambda = 1;
        while(lambda > epK)
            zdNew = zd + lambda * d;
            QNew = f(zdNew);
            if(QNew < Q)
                zd = zdNew;
                Q = QNew;
                break
            else
                lambda = lambda/2; % kontrakcja
            end
        end
        
        % --- Krok 7 ------------------------------------
        if max(abs(zdOld - zd))
            R = false;
            break
        else
            if R
                fprintf('KONIEC - KIERUNEK NAJSZYBSZEGO SPADKU I BRAK POPRAWY\n');
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
iter
zd
Q
end