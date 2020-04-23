function [ CF ] = correction( Load,RPM )
if Load<=1
   %CF = (-3*10^(-9))*RPM^2+5*10^(-5)*RPM+.7088;
   CF=-7.9*10^(-9)*RPM^2+.00005*RPM+.73;
   if Load<=.9
       CF=-8*10^(-9)*RPM^2+.00005*RPM+.73;
       %CF = CF-(1-Load)/4;
   end
end
 
end
