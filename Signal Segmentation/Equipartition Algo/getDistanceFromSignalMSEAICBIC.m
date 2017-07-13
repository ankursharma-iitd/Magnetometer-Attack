%Computes the distance matrix for EP using
%mode = 1:  FFT
%mode = 2:  Polynomial method
%mode = 3:  Wavelet method
%mode = 4:  FFT + Polynomial method (basis selection) + Wavelet method
%Sapo : min number of coefficients  for the FFT 
%Seos : max number of coefficients for the FFT 
%2S(1)-1: number of coefficients Polynomial
%CRITERION: parameter for the FFT

function [d,meth,Smeth] =  getDistanceFromSignalMSEAICBIC(Y,Sapo,Seos,CRITERION,AICBIC_)

Y = (Y-mean(Y));
Y = Y/sqrt(var(Y));

M = length(Y);
S = Sapo;
d(1:M,1:M) = 0;
meth(1:M,1:M) = 0;
Smeth(1:M,1:M) = 1;

for i=1:M,
    i/M
    for j=i+1:M,
        AIC_ = 10000000*ones(Seos,3);
        MSE = 10000000*ones(Seos,3);
        for S=Sapo:Seos,
            s = Y(i:j);
            L = length(s);
            S_real = min(S,ceil(L/2));
            [s_,errorF,errorTime] = reconstructFFTMinEnergy(s,S_real,CRITERION);
            errorF = mean((s-s_).^2)+0.0001*S_real;
            S_real = min(2*S-1,ceil(L/2));
            [s_,errorP,errorTime] = reconstructPolyMinEnergy(s',S_real);
            errorP = mean((s'-s_).^2)+0.0001*S_real/2;
            S_real = min(2*S,ceil(L/2));
            [s_,errorW,errorTime] = reconstructWaveletsMinEnergy(s',S_real);
            errorW = mean((s-s_').^2)+0.0001*S_real/2;
            
            if AICBIC_ == 1,
                AIC_(S,1) = 2*S_real+L*log(errorF);
                AIC_(S,2) = 2*S_real+L*log(errorP);
                AIC_(S,3) = 2*S_real+L*log(errorW);
            else %BIC_ log(e)+log(n)*K/n
                AIC_(S,1) = (log(L)*S_real/L) +log(errorF);
                AIC_(S,2) = (log(L)*S_real/L) +log(errorP);
                AIC_(S,3) = (log(L)*S_real/L) +log(errorW);
            end
            
            MSE(S,1) = errorF;
            MSE(S,2) = errorP;
            MSE(S,3) = errorW;
        end
        
        temp = min(min(AIC_));
        [S,method] = find(AIC_ == temp);
        S = S(1);
        method = method(1);
        
        meth(i,j) = method;
        meth(j,i) = method;
        Smeth(i,j) = S;
        Smeth(j,i) = S;
        
        d(i,j) = MSE(S,method);
        d(j,i) = MSE(S,method);
        
    end
end


figure;
imagesc(d);

figure;
imagesc(meth);

figure;
imagesc(Smeth);
