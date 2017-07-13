
%Import patient data
%Inputs are the 'filename', 'number of samples to segment', 'sample rate', 'threshold value'
patientData = importdata('sample.csv'); 
n = length(patientData);
Time = patientData(1:n,1);
EEG1 = patientData(1:n,2);
length_signal = 1166; %number of samples I want to segment
EEGSignal1_42sec = EEG1(1:length_signal);
fs = length_signal/(Time(length_signal)-Time(1)); %better sample rate, better results. You can manually enter the sample rate too.

lim = 5; %number of signals to be taken while detecting false segments
threshold = 1; %threshold for detecting false segments. Can be manually adjusted
winlength = 17; %window length to be taken

envelope(EEG1(1:length_signal),50,'rms'); %prints the envelope
[yupper, ylower] = envelope(EEG1(1:length_signal),50,'rms'); %getting the envelope for detecting segment checkers
%I have used yupper envelope in my case

% %Perform adpative signal segmentation
[~,X_m,Y_m] = ga2( yupper, Time, 1,0,length_signal, fs, winlength, EEGSignal1_42sec);
FitnessFunction = @(x) fitness(x,X_m,Y_m, length_signal);
numberOfVariables = 2;
[x,fval] = ga(FitnessFunction,numberOfVariables);
[savGolSegments,~,~] = ga2( yupper, Time, x(1),x(2),length_signal, fs, winlength, EEGSignal1_42sec); %x(1) - A1, x(2) - F1
% EEGSignal1_42sec = yupper;
%Plot originalclea EEG signal
figure; 
hold;
plot(Time(1:length_signal),EEGSignal1_42sec);
xlabel('Time(seconds)');
ylabel('Amplitude');
title('Original Signal');

time = Time(1:length_signal);
timeseg = [];
for idx = 1:size(savGolSegments,2)
    curr_seg = savGolSegments(idx);
    for idx2 = 1:size(time,1)
        if(time(idx2) == curr_seg)
            timeseg = [timeseg idx2];
        end
    end
end
% yupper = sgolayfilt(yupper, 3 , 51);
[peaks,locs] = findpeaks(yupper);
btsegments = [];
means = [];
for idx = 1:size(savGolSegments,2)
   disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------'); 
    curr_sample = timeseg(idx);
    curr_seg = savGolSegments(idx);
    disp('CURRENT SAVGOL SAMPLE');
    disp(curr_sample);
    
    selected_y = [];
    selected_x = [];
    x_left = [];
    x_right = [];
    check1 = 0;
    check2 = 0;
    for idx3 = 1:(size(locs,1)-lim)
        if((locs(idx3+lim-1)<(curr_sample)) && (locs(idx3+lim)>(curr_sample)))
            check1 = idx3;
            for idx4 = 1:lim
                x_left = [x_left peaks(check1+idx4-1)];
            end
            for idx5 = 1:lim
                x_right = [x_right peaks(check1+lim+idx5-1)];
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
                x_left = [x_left peaks(check2+idx5-1)];
                x_right = [x_right peaks(check2+lim+idx5)];               
            end
            disp('MATCHED CASE');
            disp('SAMPLE VALUES TO THE LEFT :-');
            disp(x_left);
            disp('SAMPLE VALUES TO THE RIGHT :-');
            disp(x_right);
            break;
        end       
    end

    amp_mean1 = mean(x_left);
    amp_mean2 = mean(x_right);
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

%Plot Savitzky-Golay filtered EEG signal
figure; 
hold;
plot(Time(1:length_signal),EEGSignal1_42sec);
% plot(Time(1:length_signal), EEGSignal1_42sec);
title('Normal Signal');
xlabel('Time(seconds)');
ylabel('Amplitude');
for idx = 1:size(savGolSegments,2)
plot([savGolSegments(idx) savGolSegments(idx)], [-600 -400],'g');
end
for idx = 1:size(btsegments,2)
plot([btsegments(idx) btsegments(idx)], [-600 -400],'r');
end

disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
disp(means);
disp('ORIGINAL SEGMENTS');
disp(savGolSegments);
disp('ACTUAL SEGMENTS');
disp(btsegments);
disp(x);
disp(fval);

