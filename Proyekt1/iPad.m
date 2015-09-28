% ITAM
% Análisis Aplicado con el doctor Zeferino Parada
% 14. september 2015
% Morales Mendoza, Rodrigo Andrés 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I-Pad Logistische Funktion - Wachstums Modell
%       Curva de crecimiento logístico
clc;
% clear;
% https://www.wikiwand.com/de/IPad#/Verkaufszahlen
% Datos trimestrales de ventas para IPAD:
% Kalenderjahr	1. Quartal 	2. Quartal	3. Quartal	   4. Quartal	     Summe
            % (Jan.–März)  (Apr.–Juni)   (Juli–Sep.)    (Okt.–Dez.)
% 2011   4.694.000[56]	9.246.000[57]	11.123.000[58]	15.434.000[59]	40.497.000
% 2012	11.798.000[60]	17.042.000[61]	14.036.000[62]	22.860.000[63]	65.736.000
% 2013	19.477.000[64]	14.617.000[65]	14.079.000[66]	26.035.000[67]	74.208.000
% 2014	16.350.000[68]	13.276.000[69]	12.316.000[70]	21.419.000[71]  63.361.000
% 2015	12.623.000[72]	10.931.000[73]	23.554.000
% Summe................................................................ 282.145.000
% Sólo tomamos 16 datos, como se indica en el proyecto
global v t;
v = [4694, 9246, 11123, 15434,...
    11798, 17042, 14036, 22860,...
    19477, 14617, 14079, 26035,...
    16350, 13276, 12316, 21419]';
v = v/1000; % trabajando en millones (de ventas)
n = length(v);
t = (1:n)';
% x=[.4081 19.2104 2.0526]'; solución real*
x0 = [0.01 500 3]'; % fürs I-Pad
% x0 = [0.005 1 30]'; % für Agrikultur
error0 = F(x0);
tic
%[X,iter] = busqLin('F',x0) % CORE del sistema
%[X,iter] = MetodoNewton('F',x0) %NO sirvió
%[X,iter] = MetodoNewton2('F',x0)
[X,iter] = MetodoNewton3('F',x0)
%[X,iter] = MetodoMaximoDesc('F',x0) %
toc

time = linspace(1,n,100);
fit0 = fi(x0(1),x0(2),x0(3),time);
fit2 = fi(X(1),X(2),X(3),time);
% Plots
figure(11)
plot(t,v,'*k',time,fit0,'-r')
legend('Data','1st Try')
title('I-Pad')
xlabel('Tiempo')
ylabel('Ventas en millones')
hold on;
plot(time,fit2,'-b')
legend('Data','1st Try','FitAjustado')
hold off;
