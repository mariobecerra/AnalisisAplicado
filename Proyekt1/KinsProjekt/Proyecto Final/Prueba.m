%%PRUEBA
%Ipad
x = [0.01,30,3]';

[parametrosIpad,iterIpad,fxIpad,normaIpad]=LevenbergMarquadt('resipad',x);


% %Trigo

x = [0.005,30,8]';

[parametrosTrigo,iterTrigo,fxTrigo,normaTrigo]=LevenbergMarquadt('restrigo',x);

