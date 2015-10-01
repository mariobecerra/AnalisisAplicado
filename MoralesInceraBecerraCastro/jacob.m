% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Andrés Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Curva de crecimiento logístico
% función de matriz Jacobiana
% iPad y Trigo
function [ Ja ] = jacob( fname,x)
% Encuentra la matriz Jacobiana para Levenberg-Marquardt,
% como se propone y explica en Nocedal, cap. 10.
h=1.e-6;
n=length(x);
Ja=zeros(length(feval(fname,x)),n);

for j = 1:length(x)
    e = zeros(n,1);
    e(j) = 1;
    fxh = feval(fname,x+h*e);
    fx = feval(fname,x);
    Ja(:,j) = (1/h)*(fxh-fx);
end
