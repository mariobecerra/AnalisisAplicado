% ITAM, Septiembre 2015
% An�lisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% Jos� Carlos Castro 127049
% Jos� Manuel Incera 125360
% Rodrigo Andr�s Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curva de crecimiento log�stico
% Agrikultur Logistische Funktion - Wachstums Modell
% Agricultura

clear;
global VECTOR_DE_DATOS v t;
VECTOR_DE_DATOS = [11.72 13.38 14.10 13.87 14.80 15.68 14.36 16.30 16.91 18.16 18.43 18.70 20.46 19.16 20.02 22.41 21.21 22.81 23.97 23.27 23.80 25.59 24.93 26.59]';
v = VECTOR_DE_DATOS;
n = length(VECTOR_DE_DATOS);
x0 = [0.005 1 30]';
t = (1:n)';
tic
[X, iter] = GaussNewton('residuales',  1.e-06, 10000, x0);
toc

disp('La soluci�n final para el problema de trigo es:')
disp(strcat('r=', num2str(X(1))));
disp(strcat('K=', num2str(X(2))));
disp(strcat('P0=', num2str(X(3))));

time = linspace(1,n,100);
fit = func_log(X(1),X(2),X(3),time);

% Plots
figure(11)
plot(t,VECTOR_DE_DATOS,'*k')
title('Agricultura','FontSize',18)
xlabel('Tiempo','FontSize',18)
ylabel('Producci�n de trigo por hect�rea de cultivo','FontSize',18)
hold on;
plot(time,fit,'-b')
hold off;

time3 = linspace(-100,100,1000);
fit3 = func_log(X(1),X(2),X(3),time3);
figure(13) %para ver el ss de nuestra funci�n
plot(time3,fit3,'-b')
title('Producci�n estacionaria Trigo - Predicci�n del modelo','FontSize',18)
xlabel('Tiempo','FontSize',18)
ylabel('Producci�n de trigo por hect�rea de cultivo','FontSize',18)
