% ITAM
% Análisis Aplicado con el doctor Zeferino Parada
% 14. september 2015
% Morales Mendoza, Rodrigo Andrés 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ipad Logistische Funktion
% Esta función de error es la que se desea minimizar.
function [ FR ] = F( x )
% funcion de suma de errores
% F = (1/2)* Sum_i=1^n (fi -vi)^2
%   donde fi es la logistische Funktion
% In:
%   x(1) = r; x(2)=K; x(3) =P0
%   v son los datos de trigo o Ipad
% Out:
%   fx es el error total
global v t;
FR = 1/2 * (fi(x(1),x(2),x(3),t)-v)'*(fi(x(1),x(2),x(3),t)-v);
end
