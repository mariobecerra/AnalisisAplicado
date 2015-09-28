function [ trigo ] = restrigo(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
trigo=nan(24,1);
y=nan(24,1);
y(1) = 11.72;
y(2) = 13.38;
y(3) = 14.1;
y(4) = 13.87;
y(5) = 14.8;
y(6) = 15.68;
y(7) = 14.36;
y(8) = 16.3;
y(9) = 16.91;
y(10) = 18.16;
y(11) = 18.43;
y(12) = 18.7;
y(13) = 20.46;
y(14) = 19.16;
y(15) = 20.02;
y(16) = 22.41;
y(17) = 21.21;
y(18) = 22.81;
y(19) = 23.97;
y(20) = 23.27;
y(21) = 23.8;
y(22) = 25.59;
y(23) = 24.93;
y(24) = 26.59;
for i=1:length(y)
    
  trigo(i)=x(2)/(1+   (x(2)/x(3)-1)*exp(-x(1)*(i-1)))-y(i);

end


end

