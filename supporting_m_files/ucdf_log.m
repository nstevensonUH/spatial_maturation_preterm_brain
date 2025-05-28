function [b,n] = ucdf_log(x, N);

x1 = log10(x);
th = linspace(min(x1), max(x1), N);
b = zeros(1,N);
for ii = 1:N;
    b(ii) = sum(x1>th(ii));
end
%b = b./length(x1);
n = th-(mean(diff(th))/2);
