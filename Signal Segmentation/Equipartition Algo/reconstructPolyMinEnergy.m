%polymonial fitting 
function [y,error,errorTime] =  reconstructPolyMinEnergy(s,coeff)

n = coeff;
L = length(s);
t = [0:1:L-1]/L;
p = polyfit(t,s,n);
y = polyval(p,t);




 errorTime = norm(s-y)/L;
 
 energy = sum(y.^2);
 totalEnergy = sum(s.^2);
 
 error = abs(energy-totalEnergy)/L;
 %error = errorTime;

