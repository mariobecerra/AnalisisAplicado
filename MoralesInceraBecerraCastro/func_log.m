% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Andrés Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Curva de crecimiento logístico
% función de aproximación fi
% iPad y Trigo

function [ fi ] = fi( r, k, Po, ti )
	fi = (k.*Po.*exp(r.*ti))./(Po.*exp(r.*ti)+k-Po);
end

