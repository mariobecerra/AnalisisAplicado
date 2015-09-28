% ITAM
% Análisis Aplicado con el doctor Zeferino Parada
% 14. september 2015
% Morales Mendoza, Rodrigo Andrés 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de solución 
function [ x, it ] = busqLin( fu,x0 )
global v t; %en caso de que se requieran
% In: función, x0
% Out: x, iteraciones (arg min f(x)) y num de iteraciones
it = 0;
x=x0;
g  = gradiente(fu,x);
ng = norm(g);

end

