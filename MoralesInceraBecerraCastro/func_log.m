% ITAM, Septiembre 2015
% An�lisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% Jos� Carlos Castro 127049
% Jos� Manuel Incera 125360
% Rodrigo Andr�s Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Curva de crecimiento log�stico
% funci�n de aproximaci�n fi
% iPad y Trigo

function [ fi ] = fi( r, k, Po, ti )
	fi = (k.*Po.*exp(r.*ti))./(Po.*exp(r.*ti)+k-Po);
end

