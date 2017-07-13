function y = fitness(x, X_m, Y_m, length_signal)
length_sig = size(X_m,2);
X_m = X_m(1:length_sig-1);
Y_m = Y_m(1:length_sig-1);
mean_x = mean(X_m);
mean_y = mean(Y_m);
y = 0;
for idx = 1:length_sig-1
    y = y + (abs(ceil(x(1)*(X_m(idx) - mean_x) + x(2)*(Y_m(idx) - mean_y))))^2;   
end
y = y/length_signal;
end