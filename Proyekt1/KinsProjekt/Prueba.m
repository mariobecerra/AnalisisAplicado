%%PRUEBA
%Ipad

timetot=0;
nn=1000
for i=1:nn
x = [0.01,30,3]';

tic
[parametrosIpad,iterIpad,fxIpad,normaIpad]=LevenbergMarquadt('resipad',x);


% %Trigo

x = [0.005,30,8]';

[parametrosTrigo,iterTrigo,fxTrigo,normaTrigo]=LevenbergMarquadt('restrigo',x);

time =toc
timetot=time+timetot;
end
timeprom=timetot/nn