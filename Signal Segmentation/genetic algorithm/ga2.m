%This function segments the signal on the basis of Modified Varri Method
%Please refer tho the research paper on 'modified signal segmentation'
%before doing this
function [savGolSegments X_m Y_m] = ga2(signalValue, time,A1,F1, length_signal, fs, winlength, actual_signal)
%Savitzky-Golay Filter
% sgfInput = sgolayfilt(eegSignal, 3 , 51);
% to test without filer, uncomment the line below and comment above
input_to_MV = signalValue;
input_to_MV = input_to_MV.'; %did transpose to match dimensions 
input_2_to_MV = actual_signal;
input_2_to_MV = input_2_to_MV.';

%Calculate the G variations for the filtered signal
[G_savGol X_m Y_m] = modifiedVarri(input_to_MV, A1, F1, length_signal, winlength, input_2_to_MV);

%G threshold: mean(G)
threshold_savGol = mean(G_savGol); %threshold used for G_m

%Find local maxima
[maxima_savGol, locMax_savGol]= findpeaks(G_savGol);

gScaledTime = time(1:fs:end); %scaling the time frame back
savGolSegments = [];
for idx = 1:size(maxima_savGol,2)
    if(maxima_savGol(idx) > threshold_savGol)
        savGolSegments = [savGolSegments gScaledTime(locMax_savGol(idx))];
    end
end
end

%Calculates A_dif for a given window
function[A_diff] = adif(window)
x_upper = [];
x_lower = [];
x_mid = mean(window);
for idx = 1:size(window,2)
    if(window(idx) > x_mid)
        x_upper = [x_upper window(idx)];
    else
        if(window(idx) < x_mid)
            x_lower = [x_lower window(idx)];
        end
    end
end
mean_up = mean(x_upper);
mean_down = mean(x_lower);
A_diff = abs(mean_up - mean_down);
end

%Calculates F_dif for a given window
function [F_dif] = fdif(window)
F_dif = 0;
[~, idx] = findpeaks(window);
F_dif = mean(abs(diff(idx)));
end

% %Calculates A_dif for a given window
% function[A_dif] = adif(window)
% A_dif = 0;
% for idx = 1:size(window,2)
%     A_dif = A_dif + abs(window(idx));
% end
% end
% 
% %Calculates F_dif for a given window
% function [F_dif] = fdif(window)
% F_dif = 0;
% for idx = 2:size(window, 2)
%     F_dif = F_dif + abs(window(idx) - window(idx - 1));
% end
% end

%Modified Varri Method
function [G_m, X_m, Y_m] = modifiedVarri(signal, A1, F1, lenlen, winlength, signal2)
X_m = [];
Y_m = [];
winLen = winlength; %window length (samples)
signalLen = lenlen; %total signal length (samples)
numWin = ceil(signalLen/winLen); %Num of windows for signal
winNumElem =  floor(signalLen/numWin);
G_m = zeros(1,numWin);    %initialize measure difference function
%Calculate G_m
    for win = 1:numWin-1 %handles the last window case(not necessarily same len as other windows)
    if(win == numWin-1)
        winNumElemNext = size(signal,2)/(win+1);
    else
        winNumElemNext = winNumElem;
    end
    currWin = signal((win-1)*winNumElem + 1 : win*winNumElem);
    nextWin = signal(win*winNumElem + 1: (win+1)*winNumElemNext);
    currWin2 = signal2((win-1)*winNumElem + 1 : win*winNumElem);
     nextWin2 = signal2(win*winNumElem + 1: (win+1)*winNumElemNext);
    G_m(win) = A1 * abs(adif(nextWin) - adif(currWin)) + F1 * abs(fdif(nextWin2) - fdif(currWin2));
    X_m = [X_m (abs(adif(nextWin) - adif(currWin)))];
    Y_m = [Y_m (abs(fdif(nextWin2) - fdif(currWin2)))];
    end
end