% ITAM
% An�lisis Aplicado con el profesor Zeferino Parada
% 14. august 2015
% Morales Mendoza, Rodrigo Andr�s 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ vgradiente ] = gradiente( fname, x )
%Aproximaci�n num�rica al gracdiente de fname: R^n --> R
% en x usando diferencias hacia adelante 
h = 1.e-05;
n = length(x);
fx = feval(fname,x);
vgradiente=zeros(n,1);
for k=1:n
    xaux =x;
    xaux(k) = xaux(k) + h;
    fx2 = feval(fname,xaux);
    vgradiente(k) = (fx2-fx)/h;
end

