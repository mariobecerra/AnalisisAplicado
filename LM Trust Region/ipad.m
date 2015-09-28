% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curva de crecimiento logístico
% IPAD

clear;
global VECTOR_DE_DATOS;
VECTOR_DE_DATOS = [3.27, 4.19, 7.33, 4.69, 9.25, 11.12, 15.3, 11.8, 17, 14, 22.9, 19.5, 14.6, 14.10, 26, 16.35]'; % trabajando en millones (de ventas)
n = length(VECTOR_DE_DATOS);
x0 = [0.01 500 3]';
t = (1:n)';
tic
[X, iter] = GaussNewton('residuales',  1.e-06, 10000, x0);
toc

disp('La solución final para el problema de iPads es:')
disp(strcat('r=', num2str(X(1))));
disp(strcat('K=', num2str(X(2))));
disp(strcat('P0=', num2str(X(3))));

time = linspace(1,n,100);
fit = func_log(X(1),X(2),X(3),time);

% Plots
figure(12)
plot(t,VECTOR_DE_DATOS,'*k')
title('iPad','FontSize',18)
xlabel('Tiempo','FontSize',18)
ylabel('Ventas en millones','FontSize',18)
hold on;
plot(time,fit,'-b')
hold on
hold off;
