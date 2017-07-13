%This code is a simple (not speed optimized) implementation of EP   
% based Simultaneous Segmentation and Modelling of Signals, see papers [1-3]. 
%
% [1] C. Panagiotakis, K. Athanassopoulos and G. Tziritas, The equipartition of curves,
%  Computational Geometry: Theory and Applications, Vol. 42, No. 6-7, pp. 677-689, 2009.
% 
% [2] C. Panagiotakis and G. Tziritas, Signal Segmentation and Modelling based on Equipartition Principle, 
% International Conference on Digital Signal Processing, 2009.
% 
% [3] C. Panagiotakis and G. Tziritas, Simultaneous Segmentation and Modelling of Signals based 
% on an Equipartition Principle, 20th International Conference for Pattern Recognition  (ICPR), 2010.



close all;
clear all;

S = 7;%number of coefficients 
N = 3;%number of segments 
%mode = 1;  %FFT
%mode = 2;  %Polynomial method
%mode = 3;  %Wavelet method
mode = 4;  %Selection between FFT + Polynomial method + Wavelet (basis selection)

file = 'data\sample.txt'; %data file

%read data
fid = fopen(file, 'r');
a = fscanf(fid, '%g', [inf]);    
a = a';
fclose(fid);

y = a';

figure;
plot(y);
xlabel('samples');
title('signal');
%compute distance matrix 
warning('off','all');
[d,meth] = getDistanceFromSignal(y,4,S,1);

%write distance matrix 

fileD = sprintf('%s_D.txt',file(1:length(file)-4));
fid = fopen(fileD, 'w');

for i=1:size(d,1),
    fprintf(fid,'%f ',d(i,:));
    fprintf(fid,'\n');
end
fclose(fid);


%%%EP starts

command = sprintf('ep %s %d data\\ epRes_signal4SegN1-noise',fileD,N); 
dos(command);
[Sol1,error1] =  getBestSolution(d,'solutionsEQP.txt',N,1);
[Sol2,error2] =  getBestSolution(d,'solutions.txt',N,2);
if error1 < error2,
    Sol = round(Sol1);
else
    Sol = round(Sol2);
end
fclose(fid);

%Sol stores the best EP solution 


Segments = Sol;

[y_,error,method] =  getReconstruction(y,Segments,mode,S,1);

figure;
hold on;
plot(y,'LineWidth',1.2);
hold on;
plot(y_,'r');
hold on;
%plot(x(points),0,'gs','MarkerFaceColor','g','MarkerSize',6);
for i=2:length(Segments)-1,
    plot([Segments(i) Segments(i)],[min(y) max(y)],'-.k');
end
for i=1:length(Segments)-1,
    if method(i) == 1,
        text(-10+(Segments(i)+Segments(i+1))/2,(min(y)+max(y))/2,'Fourier');
    elseif method(i) == 2,
        text(-10+(Segments(i)+Segments(i+1))/2,(min(y)+max(y))/2,'Polynomial');
    else
        text(-10+(Segments(i)+Segments(i+1))/2,(min(y)+max(y))/2,'Wavelets');
    end 
end

s = sprintf('N = %d,  Error = %2.4f, S = %d ',N,error, S);
title(s); 
legend('Original','Reconstruction','Segmentation');
xlabel('samples');




