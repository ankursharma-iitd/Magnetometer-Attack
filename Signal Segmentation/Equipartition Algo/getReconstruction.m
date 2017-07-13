%Computes the distance matrix for EP using
%mode = 1:  FFT 
%mode = 2:  Polynomial method
%mode = 3:  FFT + Polynomial method (basis selection)
%segments tmhmatopoihsh
%S : # suntelestwn FFT
%2S-1: # suntelestwn Polynomial
%CRITERION: parametros gia ton FFT

function [z,MaxError,Method] =  getReconstruction(Y,segments,mode,S,CRITERION)

M = length(Y);
N = length(segments);
MaxError = 0;
d(1:M,1:M) = 0;
meth(1:M,1:M) = 0;
z = [];
Method = [];
if mode == 1,
    for i=1:N-1,
        apo = segments(i);
        eos = segments(i+1);
        if i == N-1,
            eos = eos+1;
        end
        s = Y(apo:eos-1);
        L = length(s);
        S_real = min(S,ceil(L/2));
        [s_,error,errorTime,mf] = reconstructFFTMinEnergy(s,S_real,CRITERION);
        z = [z s_];
        MaxError = max(MaxError,error);
    end
elseif mode == 2,
    for i=1:N-1,
        apo = segments(i);
        eos = segments(i+1);
        if i == N-1,
            eos = eos+1;
        end
        s = Y(apo:eos-1);
        L = length(s);
        S_real = min(2*S-2,ceil(L/2));
        [s_,error,errorTime] = reconstructPolyMinEnergy(s',S_real);
        MaxError = max(MaxError,error);
        z = [z s_];
    end
elseif mode == 3,
    for i=1:N-1,
        apo = segments(i);
        eos = segments(i+1);
        if i == N-1,
            eos = eos+1;
        end
        s = Y(apo:eos-1);
        L = length(s);
        S_real = min(2*S-1,ceil(L/2));
        [s_,error,errorTime] = reconstructWaveletsMinEnergy(s',S_real);
        MaxError = max(MaxError,error);
        z = [z s_];
    end
else
    for i=1:N-1,
        apo = segments(i);
        eos = segments(i+1);
        if i == N-1,
            eos = eos+1;
        end
        s = Y(apo:eos-1);
        L = length(s);
        
        S_real = min(S,ceil(L/2));
        [s_F,errorF,errorTime] = reconstructFFTMinEnergy(s,S_real,CRITERION);
        s_F = s_F';
        errorF = errorF+0.0001*S_real;
        S_real = min(2*S-1,ceil(L/2));
        [s_P,errorP,errorTime] = reconstructPolyMinEnergy(s',S_real);
        S_real = min(2*S,ceil(L/2));
        errorP = errorP+0.0001*S_real/2;
        [s_W,errorW,errorTime] = reconstructWaveletsMinEnergy(s',S_real); 
        errorW = errorW+0.0001*S_real/2;
        
        [error method] = min([errorF errorP errorW]);
        
        MaxError = max(MaxError,error);
        s_ = s_F;
        if method == 2,
            s_ = s_P;
        elseif method == 3,
            s_ = s_W;
        end
        Method = [Method method];
        z = [z s_];
    end
 end



