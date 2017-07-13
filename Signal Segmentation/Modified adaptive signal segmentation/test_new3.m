%% My implementation of 
%
%An Improved Signal Segmentation Using Moving Average and Savitzky-Golay Filter
%Azami et al. 2012

function [freqSegments] = test_new3(eegSignal, time,A1,F1, length_signal, fs, winlength)
%Savitzky-Golay Filter
% sgfInput = sgolayfilt(eegSignal, 3 , 5);
% to test without filer, uncomment the line below and comment above
sgfInput = eegSignal;
sgfInput = sgfInput.'; %did transpose to match dimensions of tsmovavg 

%Calculate the G variations for the filtered signal
G_savGol = modifiedVarri(sgfInput, A1, F1, length_signal, winlength);
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
freqSegments = [];
for idx = 1:size(maxima_savGol,2)
    if(maxima_savGol(idx) > threshold_savGol)
        freqSegments = [freqSegments gScaledTime(locMax_savGol(idx))];
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
% disp(A_diff);
end

%Calculates F_dif for a given window
function [F_dif] = fdif(window)
F_dif = 0;
[~, idx] = findpeaks(window);
F_dif = mean(diff(idx));
end
%Modified Varri Method
function [G_m] = modifiedVarri(signal, A1, F1, lenlen, winlength)
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
      disp('blehhhhh start');
      disp(adif(nextWin));
      disp(adif(currWin));
      if(F1 == 0) 
           btt = (A1*(abs(adif(nextWin) - adif(currWin))));
      else
          if(A1 == 0)
             btt =(F1*(abs(fdif(nextWin) - fdif(currWin)))); 
          else
              x = abs(adif(nextWin) - adif(currWin));
              y = abs(fdif(nextWin) - fdif(currWin));
              disp('X down');
              disp(x);
              disp('Y down');
              disp(y);
              btt = (A1*x) + (F1*y);
          end
      end
      G_m(win) = btt;
      disp(win);
      disp(btt);
      disp(G_m(win));
      disp('blehhhhh end');
end
end