% ITAM
% An�lisis Aplicado con el doctor Zeferino Parada
% 14. september 2015
% Morales Mendoza, Rodrigo Andr�s 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I-Pad Logistische Funktion - Wachstums Modell
%       Curva de crecimiento log�stico
% AGRICULTURA
clc;
% clear;
global v t;
v = [11.72 13.38 14.10 13.87 14.80 15.68 14.36 16.30 16.91 18.16 18.43 18.70 20.46...
19.16 20.02 22.41 21.21 22.81 23.97 23.27 23.80 25.59 24.93 26.59]';
v = v/1000; % trabajando en millones
n = length(v);
t = (1:n)';
% x=[.4081 19.2104 2.0526]'; soluci�n real
% x0 = [0.01 500 3]'; % f�rs I-Pad
x0 = [0.005 1 30]'; % f�r Agrikultur
error0 = F(x0);
tic
%[X,iter] = busqLin('F',x0) % CORE del sistema
%[X,iter] = MetodoNewton('F',x0) %NO sirvi�
%[X,iter] = MetodoMaximoDesc('F',x0) %
[X,iter] = Metodo_Descenso_Newton('F',x0)
toc

time = linspace(1,n,100);
fit0 = fi(x0(1),x0(2),x0(3),time);
fit2 = fi(X(1),X(2),X(3),time);
% Plots
figure(12)
%plot(t,v,'*k')
plot(t,v,'*k',time,fit0,'-r')
legend('Data','1st Try')
title('I-Pad')
xlabel('Tiempo')
ylabel('Ventas en millones')
hold on;
plot(time,fit2,'-b')
legend('Data','1st Try','FitAjustado')
hold off;
