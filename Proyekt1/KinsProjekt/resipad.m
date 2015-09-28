function [ ipad ] = resipad(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
ipad=nan(16,1);
y = [3.27, 4.19, 7.33, 4.69, 9.25, 11.12, 15.3, 11.8, 17, 14, 22.9, 19.5, 14.6, 14.10, 26, 16.35]';

for i=1:length(y)
    
  ipad(i)=x(2)/(1+   (x(2)/x(3)-1)*exp(-x(1)*(i-1)))-y(i);

end


end

