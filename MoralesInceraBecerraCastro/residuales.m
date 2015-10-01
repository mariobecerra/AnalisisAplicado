% ITAM, Septiembre 2015
% An�lisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% Jos� Carlos Castro 127049
% Jos� Manuel Incera 125360
% Rodrigo Andr�s Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Curva de crecimiento log�stico
% funci�n de residuales
% iPad y Trigo

function [ res ] = residuales(x)
% Encuentra los Residuales entre la aproximacion y las ventas reales
% Par�metros:
% 	x = vector de tres par�metros donde:
% 		x(1) = r, tasa de aceptaci�n
%		x(2) = K, cantidad m�x permitida en la poblaci�n
%		x(3) = P0, poblaci�n al tiempo 0
%	v = variable global, vector con ventas de iPad o prod trigo

	global VECTOR_DE_DATOS;

	u = 0:(length(VECTOR_DE_DATOS) - 1);
    u = u';
	res = x(2) ./ ( 1 + (x(2) / x(3) - 1) .* exp( - x(1) .* u )) - VECTOR_DE_DATOS;

end

