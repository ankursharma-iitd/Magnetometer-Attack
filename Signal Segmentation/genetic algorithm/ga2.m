%% My implementation of 
%
%An Improved Signal Segmentation Using Moving Average and Savitzky-Golay Filter
%Azami et al. 2012

function [savGolSegments X_m Y_m] = ga2(eegSignal, time,A1,F1, length_signal, fs, winlength, EEGSignal1_42sec)
%Savitzky-Golay Filter
% sgfInput = sgolayfilt(eegSignal, 3 , 51);
% to test without filer, uncomment the line below and comment above
sgfInput = eegSignal;
sgfInput = sgfInput.'; %did transpose to match dimensions of tsmovavg 
sgfInput2 = EEGSignal1_42sec;
sgfInput2 = sgfInput2.';

%Calculate the G variations for the filtered signal
[G_savGol X_m Y_m] = modifiedVarri(sgfInput, A1, F1, length_signal, winlength, sgfInput2);
figure;
hold on;
plot(G_savGol);
hold off;

%G threshold: mean(G)
threshold_savGol = mean(G_savGol);
disp('PLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS');
disp(G_savGol);

%Find local maxima
[maxima_savGol, locMax_savGol]= findpeaks(G_savGol);

gScaledTime = time(1:fs:end);
savGolSegments = [];
for idx = 1:size(maxima_savGol,2)
    if(maxima_savGol(idx) > threshold_savGol)
        savGolSegments = [savGolSegments gScaledTime(locMax_savGol(idx))];
    end
end
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
% disp(A_diff);
end

%Calculates F_dif for a given window
function [F_dif] = fdif(window)
F_dif = 0;
[~, idx] = findpeaks(window);
F_dif = mean(abs(diff(idx)));
end

%Modified Varri Method
function [G_m, X_m, Y_m] = modifiedVarri(signal, A1, F1, lenlen, winlength, signal2)
X_m = [];
Y_m = [];
winLen = winlength; %window length (samples)
signalLen = lenlen; %total signal length (samples)
numWin = ceil(signalLen/winLen); %Num of windows for signal
winNumElem =  floor(signalLen/numWin);
lastWinNumElem = size(signal,2) - ((numWin-1)*winNumElem); %Num elements in the last window

%A1 = 10; F1 = 7; %Modified Varri constants
G_m = zeros(1,numWin);    %initialize measure difference function

%Calculate G_m
for win = 1:numWin-1

    %handles the last window case(not necessarily same len as other windows)
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


