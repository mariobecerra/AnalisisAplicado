% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Andrés Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Curva de crecimiento logístico
% función de residuales
% iPad y Trigo

function [ res ] = residuales(x)
% Encuentra los Residuales entre la aproximacion y las ventas reales
% Parámetros:
% 	x = vector de tres parámetros donde:
% 		x(1) = r, tasa de aceptación
%		x(2) = K, cantidad máx permitida en la población
%		x(3) = P0, población al tiempo 0
%	v = variable global, vector con ventas de iPad o prod trigo

	global VECTOR_DE_DATOS;

	u = 0:(length(VECTOR_DE_DATOS) - 1);
    u = u';
	res = x(2) ./ ( 1 + (x(2) / x(3) - 1) .* exp( - x(1) .* u )) - VECTOR_DE_DATOS;

end

