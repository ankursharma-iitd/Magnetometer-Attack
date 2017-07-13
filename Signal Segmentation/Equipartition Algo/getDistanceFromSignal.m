%Computes the distance matrix d and the correspoding method matrix meth for EP using
%mode = 1:  FFT
%mode = 2:  Polynomial method
%mode = 3:  Wavelet method
%mode = 4:  FFT + Polynomial method + Wavelet (basis selection)
%S : # coefficients FFT
%2S-1: # coefficients Polynomial
%CRITERION: parameter for the FFT

function [d,meth] =  getDistanceFromSignal(Y,mode,S,CRITERION)

M = length(Y);

d(1:M,1:M) = 0;
meth(1:M,1:M) = 0;

if mode == 1,
    
    for i=1:M,
        for j=i+1:M,
            s = Y(i:j);
            L = length(s);
            S_real = min(S,ceil(L/2));
            [s_,error,errorTime,mf] = reconstructFFTMinEnergy(s,S_real,CRITERION);
            d(i,j) = error+0.001*S_real;
            d(j,i) = error+0.001*S_real;
        end
    end
elseif mode == 2,
    
    for i=1:M,
        for j=i+1:M,
            s = Y(i:j);
            L = length(s);
            S_real = min(2*S-2,ceil(L/2));
            [s_,error,errorTime] = reconstructPolyMinEnergy(s',S_real);
            d(i,j) = error+0.001*S_real/2;
            d(j,i) = error+0.001*S_real/2;
        end
    end
elseif mode == 3,
    
    for i=1:M,
        for j=i+1:M,
            s = Y(i:j);
            L = length(s);
            S_real = min(2*S,ceil(L/2));
            [s_,error,errorTime] = reconstructWaveletsMinEnergy(s',S_real);
            d(i,j) = error+0.001*S_real/2;
            d(j,i) = error+0.001*S_real/2;
        end
    end
else
    for i=1:M,
        for j=i+1:M,
            s = Y(i:j);
            L = length(s);
            S_real = min(S,ceil(L/2));
            [s_,errorF,errorTime] = reconstructFFTMinEnergy(s,S_real,CRITERION);
            errorF = errorF+0.0001*S_real;
            S_real = min(2*S-1,ceil(L/2));
            [s_,errorP,errorTime] = reconstructPolyMinEnergy(s',S_real);
            errorP = errorP+0.0001*S_real/2;
            S_real = min(2*S,ceil(L/2));
            [s_,errorW,errorTime] = reconstructWaveletsMinEnergy(s',S_real);
            errorW = errorW+0.0001*S_real/2;
            
            [error method] = min([errorF errorP errorW]);
            meth(i,j) = method;
            meth(j,i) = method;
            
            d(i,j) = error;
            d(j,i) = error;
        end
    end
end



figure;
imagesc(d)
title('Distance Matrix');
figure;
imagesc(meth)
title('Method Matrix');
