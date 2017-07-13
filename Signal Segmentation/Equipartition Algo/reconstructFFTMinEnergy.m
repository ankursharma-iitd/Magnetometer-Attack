%FFT fitting 

function [y,error,errorTime,mfreq,vfreq] =  reconstructFFTMinEnergy(s,coeff,criterion)

mfreq = 0;
vfreq = 0;
if criterion == 0,
    [y,error] =  reconstructFFT(s,coeff);
    errorTime = error;
    return;
end

S = coeff;
L = length(s);
f = fft(s);
L = length(s);
v = 0*[1:L];

if S >=  ceil(L/2)
    y = s;
    error = 0;
    errorTime = 0;
    return;
end

S_real = min(S,ceil(L/2));

fa = abs(f);
totalEnergy = sum(abs(fa.^2)) / L;

[temp,indexes] = sort(fa(1:ceil(L/2)),'descend');
mfreq = mean((indexes(1:S_real)-1).*abs(temp(1:S_real))) / (ceil(L/2)*mean(abs(temp(1:S_real))));

vfreq = mean(((((indexes(1:S_real)-1)/ceil(L/2))-mfreq).^2).*abs(temp(1:S_real))) / mean(abs(temp(1:S_real)));

v(indexes(1:S_real)) = 1;
v(L+2-indexes(1:S_real)) = 1;

[temp, pos] = min(indexes);
if temp == 1
    indexes(pos) = 1;
    v = v(1:L);
end

v = v';
fs = v.*f;
fas = abs(fs);
energy = sum(abs(fas.^2)) / L;

y = ifft(fs);
% errorTime = norm(s-y)/ norm(s);
% error = abs(energy-totalEnergy)/ totalEnergy;

 errorTime = norm(s-y)/L;
 error = abs(energy-totalEnergy)/L; %%%%???????????
 %error = abs(energy-totalEnergy);
 %error = errorTime;

