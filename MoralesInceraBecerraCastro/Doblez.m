% ITAM, Septiembre 2015
% An�lisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% Jos� Carlos Castro 127049
% Jos� Manuel Incera 125360
% Rodrigo Andr�s Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Curva de crecimiento log�stico
% funci�n de Doblez, como programada en clase
% iPad y Trigo
function [ p ] = Doblez( B, grad, delta )
% Se obtiene una p-direccn al prob Trust Region por el m�todo Doblez
% es decir, dentro de una bola de tama�o delta, la mejor direccn
    pN = B\-grad; % Newton
    pC = -grad * (grad' * grad) / (grad' * B * grad); % Cauchy
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

