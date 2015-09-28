% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Morales Mendoza 124341

function [ res ] = residuales(x)
% Residuales de la función logística y las ventas
% Parámetros:
% 	x = vector de tres parámetros donde:
% 		x(1) = r, tasa de aceptación
%		x(2) = K, constante de cantidad máxima permitida en la población
%		x(3) = P0, población al tiempo 0
%	v = variable global, es un vector con las ventas de iPad o la producción de trigo

	global VECTOR_DE_DATOS;

	u = 0:(length(VECTOR_DE_DATOS) - 1);
    u = u';
	res = x(2) ./ ( 1 + (x(2) / x(3) - 1) .* exp( - x(1) .* u )) - VECTOR_DE_DATOS;

end

