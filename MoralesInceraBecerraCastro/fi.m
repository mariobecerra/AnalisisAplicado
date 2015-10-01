% ITAM
% Análisis Aplicado con el profesor Zeferino Parada
% 14. september 2015
% Morales Mendoza, Rodrigo Andrés 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ipad Logistische Funktion
function [ fi ] = fi( r, k, Po, ti )
% fi = (k*Po*exp(r*ti))/(Po*exp(r*ti)+k-Po);
fi = (k.*Po.*exp(r.*ti))./(Po.*exp(r.*ti)+k-Po);
end

