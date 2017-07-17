%This code runs the genetic algorithm over a number of iterations, and then
%averages out the result which are my segments
total = [];

%following are the inputs that need to be provided
filename = 'sample.csv'; %specify the filename
length_signal = 766; %specify the length of the signal
limit = 5; %Sample Values to be considered while segment checking as if it is a true segment or false segment
threshold = 1; %difference in the mean amplitude values on the left and the right side of every segment, if greater than this threshold will imply a true segment
window_length = 17; %length of the window under consideration
number_of_iterations = 10; %number of times the algorithm runs

%runs the algorithm over a number of times, and stores all the values in a
%large vector 'total'
for idx =1:number_of_iterations
[values] = ga1(filename, length_signal, limit, threshold, window_length);
total = [total values];
end

%sort the total array
total = sort(total);

%below code merges the segments which are quite closeby
%input is the small window length under which the segments should merge
merge_length = 2.5; %this is the merge length. all the answer values under this local window merges, and gives the mean of the values

actualSegments = []; % 'actualSegments' denote the actual segments detected
idx = 1;
flag = 0;
umm = [];
while (idx < size(total,2))
    currVal = total(idx);
    nextVal = total(idx+1);
    if((nextVal - currVal) < merge_length)
        umm = [umm currVal];
        flag = 1;
    else
        if(flag == 1)
            flag = 0;
            umm = [umm currVal];
            actualSegments = [actualSegments mean(umm)];
            umm = [];
        else
            actualSegments = [actualSegments currVal];
        end
    end  
    idx = idx + 1;
end
if(flag == 1)
    umm = [umm nextVal];
    actualSegments = [actualSegments mean(umm)];
else
    actualSegments = [actualSegments nextVal];
end

sampleData = importdata(filename); 
n = length(sampleData);
Time = sampleData(1:n,1);
sampleValues = sampleData(1:n,2);
sampleSignal = sampleValues(1:length_signal);

%Plotting the Segmented Signal
figure; 
hold;
plot(Time(1:length_signal),sampleSignal);
title('Segmented Signal');
xlabel('Time(seconds)');
ylabel('Amplitude');

for idx = 1:size(actualSegments,2)
plot([actualSegments(idx) actualSegments(idx)], [-700 -300],'r'); %correct segments detected
end
