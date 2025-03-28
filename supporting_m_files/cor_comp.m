function out = cor_comp(r1, r2, r12, N, alpha)
%
% r1 is the higher correlation with FBA
% r2 is the lower correlation with FBA
% r12 is the correlation between axes (e.g SA rank with FO distance)
%
% test for overlapping correlations on dependent groups
% according to Zou GY. Toward Using Confidence Intervals to Compare Correlations. Psychol Methods. 2007;12:
% 399–413. doi: 10.1037/1082-989X.12.4.399
%

x = linspace(0, 3, 100000);
y = normpdf(x,0,1);
zalpha = x(find(cumsum(y).*(x(2)-x(1))>=0.5-alpha/2, 1));

% TEST CASE
%r1 = 0.396; r2 = .179; r12 = 0.088; N = 66;
% give [L,U] = [-0.093 0.517]

Z1 = 0.5*log((1+r1)./(1-r1));
Z2 = 0.5*log((1+r2)./(1-r2));

l1dash = Z1 - zalpha*sqrt(1/(N-3));
l2dash = Z2 - zalpha*sqrt(1/(N-3));

u1dash = Z1 + zalpha*sqrt(1/(N-3));
u2dash = Z2 + zalpha*sqrt(1/(N-3));

l1 = (exp(2*l1dash)-1) / (exp(2*l1dash)+1);
l2 = (exp(2*l2dash)-1) / (exp(2*l2dash)+1);

u1 =  (exp(2*u1dash)-1) / (exp(2*u1dash)+1);
u2 = (exp(2*u2dash)-1) / (exp(2*u2dash)+1);

pr1r2 = ((r12-0.5*r1*r2)*(1-r1^2-r2^2-r12^2)+r12^3) / ((1-r1^2)*(1-r2^2));
L = (r1-r2)-sqrt((r1-l1)^2+(u2-r2)^2 - 2*pr1r2*(r1-l1)*(u2-r2));
U = (r1-r2)+sqrt((u1-r1)^2+(r2-l2)^2 - 2*pr1r2*(u1-r1)*(r2-l2));

out = [L U];

