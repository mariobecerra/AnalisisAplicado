function [ p ] = Doblez( B, grad, delta )
% Se obtiene una aproximación al problema de región de confianza por el metodo doblez
%
% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Morales Mendoza 124341

    pN = B\-grad; % Dirección de Newton
    pC = -grad * (grad' * grad) / (grad' * B * grad); % Dirección de Cauchy

    pC_norm = norm(pC); pN_norm = norm(pN);
     
     if (pN_norm <= delta)
          p = pN;
     else
         if (pC_norm >= delta)
             p = delta * (pC / pC_norm);
         else
             dif = pN - pC;
             a = norm(dif)^2;
             b = 2 * pC' * dif  ;
             c = pC_norm^2 - delta^2  ;
             t = roots([a b c]);
             if (t(1) <= 0)
                 alpha = t(2);
             else
                 alpha = t(1);
             end
             p = pC + alpha * dif;
         end
     end
 
end

