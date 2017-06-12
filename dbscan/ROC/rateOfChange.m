print = csvread('side24-0_mag.csv');
print_t = print(:,1); % you can specify the range of the segment over here
print_x = print(:,2); % you can specify the range of the segment over here. Same as above

%can even use the filtered signal
%filt_x = sgloayfilt(print_x,3,5);

%find the peaks and calculate mean distance between indices
[peaks,idx] = findpeaks(print_x);
dist = mean(diff(idx));

roc = 1/dist;
disp("ROC = " + roc);