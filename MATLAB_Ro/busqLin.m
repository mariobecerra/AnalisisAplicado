% ITAM
% An�lisis Aplicado con el doctor Zeferino Parada
% 14. september 2015
% Morales Mendoza, Rodrigo Andr�s 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M�todo de soluci�n 
function [ x, it ] = busqLin( fu,x0 )
global v t; %en caso de que se requieran
% In: funci�n, x0
% Out: x, iteraciones (arg min f(x)) y num de iteraciones
it = 0;
x=x0;
g  = gradiente(fu,x);
ng = norm(g);

end

