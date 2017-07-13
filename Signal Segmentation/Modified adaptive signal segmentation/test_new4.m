%% This program plots the signal segments computed for the sample EEG data


%Import patient data
patientData = importdata('sample.csv');
n = length(patientData);
Time = patientData(1:n,1);
EEG1 = patientData(1:n,2);
length_signal = 1166; %length of the signal
EEGSignal1_42sec = EEG1(1:length_signal);
fs = length_signal/(Time(length_signal)-Time(1));
combine_seg_window_length = 20; %samples

lim = 7; %number of signlas to be taken while detecting false segments
threshold = 0.49; %threshold for detecting false segments
winlength = 17; %window length to be taken

[freqSegments] = test_new3(EEGSignal1_42sec,Time, 0 , 1 , length_signal, fs, winlength);
% EEGSignal1_42sec = yupper;
%Plot originalclea EEG signal
figure; 
hold;
plot(Time(1:length_signal),EEGSignal1_42sec);
xlabel('Time(seconds)');
ylabel('Amplitude');
title('Original EEG Signal');

time = Time(1:length_signal);
timeseg = [];
for idx = 1:size(freqSegments,2)
    curr_seg = freqSegments(idx);
    for idx2 = 1:size(time,1)
        if(time(idx2) == curr_seg)
            timeseg = [timeseg idx2];
        end
    end
end
disp(timeseg);

[peaks,locs] = findpeaks(EEGSignal1_42sec);

btsegments = [];
means = [];
for idx = 1:size(freqSegments,2)
   disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------'); 
    curr_sample = timeseg(idx);
    curr_seg = freqSegments(idx);
    disp('CURRENT SAVGOL SAMPLE');
    disp(curr_sample);
    
    x_left = [];
    x_right = [];
    check1 = 0;
    check2 = 0;
    for idx3 = 1:(size(locs,1)-lim)
        if((locs(idx3+lim-1)<(curr_sample)) && (locs(idx3+lim)>(curr_sample)))
            check1 = idx3;
            for idx4 = 1:lim
                x_left = [x_left locs(check1+idx4-1)];
            end
            for idx5 = 1:lim
                x_right = [x_right locs(check1+lim+idx5-1)];
            end
            disp('BETWEEN CASE');
            disp('IMMEDIATE LEFT SAMPLE');
            disp(locs(check1+lim-1));
            disp('IMMEDIATE RIGHT SAMPLE');
            disp(locs(check1+lim));
            disp('SAMPLE VALUES TO THE LEFT :-');
            disp(x_left);
            disp('SAMPLE VALUES TO THE RIGHT :-');
            disp(x_right);
            break;
        end
    end

    
    for idx3 = 1:(size(locs,1)-lim)
        if((locs(idx3+lim)==(curr_sample)))
            check2 = idx3;
            for idx5 = 1:lim
                x_left = [x_left locs(check2+idx5-1)];
                x_right = [x_right locs(check2+lim+idx5)];               
            end
            disp('MATCHED CASE');
            disp('SAMPLE VALUES TO THE LEFT :-');
            disp(x_left);
            disp('SAMPLE VALUES TO THE RIGHT :-');
            disp(x_right);
            break;
        end       
    end

    amp_mean1 = mean(diff(x_left));
    amp_mean2 = mean(diff(x_right));
    disp('MEAN OF LEFT SAMPLES');
    disp(amp_mean1);
    disp('MEAN OF RIGHT SAMPLES');
    disp(amp_mean2);
    disp('MEAN DIFFERENCE');
    disp(abs(amp_mean1 - amp_mean2));
    means = [means abs(amp_mean1 - amp_mean2)];
    
    if(abs(amp_mean1 - amp_mean2) > threshold)
        btsegments = [btsegments curr_seg];
        disp('***************************** TRUE SEGMENT ********************************************************');
    end
    
     if(abs(amp_mean1 - amp_mean2) <= threshold)
        disp('***************************** FALSE SEGMENT *******************************************************');
    end
    
end

figure; 
hold;
plot(Time(1:length_signal),EEGSignal1_42sec);
% plot(Time(1:length_signal), EEGSignal1_42sec);
title('Normal Signal');
xlabel('Time(seconds)');
ylabel('Amplitude');
for idx = 1:size(freqSegments,2)
plot([freqSegments(idx) freqSegments(idx)], [-600 -400],'g');
end
for idx = 1:size(btsegments,2)
plot([btsegments(idx) btsegments(idx)], [-600 -400],'r');
end
disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
disp(means);
disp('ORIGINAL SEGMENTS');
disp(freqSegments);
disp('ACTUAL SEGMENTS');
disp(btsegments);

