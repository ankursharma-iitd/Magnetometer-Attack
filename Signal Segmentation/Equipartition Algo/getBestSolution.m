%gets the best EP solution 
function [Sol,error] =  getBestSolution(d,filename,N,ignoreChars)

fid = fopen(filename, 'r');
a = fscanf(fid, '%f');
a = a(1+ignoreChars:length(a));

NumSol = length(a)/(N+1);

SolArray = zeros(NumSol,N+1);
k = 1;
for i=1:NumSol,
    for j=1:N+1
        SolArray(i,j) = a(k);
        k = k+1;
    end
end

Error = zeros(NumSol,1);

for i=1:NumSol,
    for j=1:N
        u = round(SolArray(i,j)+1);
        v = round(SolArray(i,j+1)+1);
        Error(i) = max(Error(i),d(u,v));
    end
end

[error pos] = min(Error);
Sol = SolArray(pos,:)+1;



