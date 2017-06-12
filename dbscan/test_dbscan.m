print = csvread('side24-0_mag.csv'); %Enter the correct name of the file
print_t = print(:,1);
print_x = print(:,2); 

%get the peaksans the indices at which they occur
[peaks, idx]=findpeaks(print_x)

%Can be used to specify custom indices. 
%The ones separated by comma are concatenated
new_idx = vertcat(idx(13:39),idx(59:118),idx(136:227));

%This code only has one iteration of DBSCAN
[index,isnoise] = dbscan(new_idx,0.32,2);

%write the results to a file
fil = fopen('dbscan.txt','w');
for i = 1:length(index)
    fprintf(fil,"%f, %f\n",index(i),isnoise(i)); 
end
fclose(fil);
