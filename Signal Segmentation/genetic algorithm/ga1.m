%This function runs the geentic algo, gets the best possible values for A1,and F1.
%Then, it chops off the false segments, and returns the remaining segments
%back to the 'ga_main.m'

function[values] = ga1(filename, length_signal, lim, threshold,winlength)

signalData = importdata(filename); 
n = length(signalData);
Time = signalData(1:n,1);
signal_total = signalData(1:n,2);
signalValues = signal_total(1:length_signal);
fs = length_signal/(Time(length_signal)-Time(1)); %better sample rate, better results. You can manually enter the sample rate too.

% envelope(EEG1(1:length_signal),50,'rms'); %prints the envelope
[yupper, ylower] = envelope(signal_total(1:length_signal),50,'rms'); %getting the envelope for detecting segment checkers
%I have used yupper envelope in my case

% %Perform adpative signal segmentation, and calling the genetic algorithm
% for two unknown variables
[~,X_m,Y_m] = ga2(yupper, Time, 1,0,length_signal, fs, winlength, signalValues); %get the X_m, and Y_m values
FitnessFunction = @(x) fitness(x,X_m,Y_m, length_signal); %call the fitness function using these parameters
numberOfVariables = 2;
% rng default
[x,fval,exitFlag, Output] = ga(FitnessFunction,numberOfVariables); %'x' vector gives me the best possible values of A1, and F1 used for the modified varri approach


% 'savGolSegments' denote the segments resulted after running the M.V.
% method with A1, and F1 calcualted from G.A.
% MV - Modified Varri
% GA - Genetic ALgo
[savGolSegments,~,~] = ga2(yupper, Time, x(1),x(2),length_signal, fs, winlength, signalValues); %x(1) - A1, x(2) - F1


%Segment Checker code
%For checking false segments, I have calculated the mean amplitude of some
%immediate left peaks, and some immediate right peaks
%If the mean difference betweeen these left and right amplitudes is less
%than a certain threshold (say 1), then the segment is a false segment, and
%can be deleted from the original list of segments
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

[peaks,locs] = findpeaks(yupper); %finding the peaks to be used later in my segment checker
btsegments = [];
means = [];
for idx = 1:size(savGolSegments,2)
   disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------'); 
    curr_sample = timeseg(idx);
    curr_seg = savGolSegments(idx);
    disp('CURRENT SAVGOL SAMPLE');
    disp(curr_sample);
    
    x_left = [];
    x_right = [];
    for idx3 = 1:(size(locs,1)-lim)
        if((locs(idx3+lim-1)<(curr_sample)) && (locs(idx3+lim)>(curr_sample))) %case 1
            check1 = idx3;
            for idx4 = 1:lim
                x_left = [x_left peaks(check1+idx4-1)]; %get the left peaks
            end
            for idx5 = 1:lim
                x_right = [x_right peaks(check1+lim+idx5-1)]; %get the right peaks
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
        if((locs(idx3+lim)==(curr_sample))) %case 2
            check2 = idx3;
            for idx5 = 1:lim
                x_left = [x_left peaks(check2+idx5-1)]; %get the left peaks
                x_right = [x_right peaks(check2+lim+idx5)]; %get the right peaks            
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
    disp(abs(amp_mean1 - amp_mean2)); %get the absolute difference between the left and right amps
    means = [means abs(amp_mean1 - amp_mean2)];
    
    if(abs(amp_mean1 - amp_mean2) > threshold)
        btsegments = [btsegments curr_seg];
        disp('***************************** TRUE SEGMENT ********************************************************');
    end
    
     if(abs(amp_mean1 - amp_mean2) <= threshold)
        disp('***************************** FALSE SEGMENT *******************************************************');
    end
    
end

%Such a segment checker can also be designed for frequency

disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
disp(means);
disp('ORIGINAL SEGMENTS');
disp(savGolSegments);
disp('ACTUAL SEGMENTS');
disp(btsegments);
disp(x); %Values of A1, and F1 at which I get the best results
disp(Output.generations); %the number of generations
disp(Output.funccount); %function counts
values = btsegments; %return the values
end

