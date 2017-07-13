%Wavelets fitting 


function [y,error,errorTime] =  reconstructWaveletsMinEnergy(s,coeff)


wname = 'db1';
n = coeff;
Ls = length(s);
[C,L] = wavedec(s,3,wname);


% weights(1:L(1)) = 1;
% apo = L(1);
% for i=2:length(L)-1,
%  weights(apo+1:apo+L(i)) = 0.5/(2^(i-2));
%  apo = apo+L(i);
% end

%[a,b] = sort(abs(C.*weights),'descend');

[a,b] = sort(abs(C),'descend');

for i=coeff+1:length(C),
    C(b(i)) = 0;
end



y = waverec(C,L,wname);

errorTime = norm(s-y)/Ls;
 
energy = sum(y.^2);
totalEnergy = sum(s.^2);
 
error = abs(energy-totalEnergy)/Ls;
%error = errorTime;

