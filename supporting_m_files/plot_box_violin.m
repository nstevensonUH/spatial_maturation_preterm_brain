function h = plot_box_violin(x, y, c, val, sd)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%figure;
hold on;
lim1 = quantile(y, [0.25 0.5 0.75]);
lim2 = [min(y(find(y>=lim1(1)-1.5*(lim1(3)-lim1(1))))) max(y(find(y<=lim1(3)+1.5*(lim1(3)-lim1(1)))))];
ols = [y(find(y<lim1(1)-1.5*(lim1(3)-lim1(1))))' y(find(y>lim1(3)+1.5*(lim1(3)-lim1(1))))'];
lim = [lim2(1) lim1 lim2(2)];

if val(2)==1
    dc = 0.25+c; if max(dc)>1; dc = dc./max(dc); end
    plot(x+rand(length(y),1)*0.25*sd+sd*0.3, y, '.', 'color', dc);
    xx = linspace(min(y), max(y), 1000);
    z = ksdensity(y, xx); z = z./max(z)*0.5;
    plot(sd*z+x+sd*0.3, xx, 'color', c, 'linewidth', 2);
end

h = zeros(length(lim), 1);
for ii = 1:length(lim)
h(ii) = plot([x-0.25 x+0.25], [lim(ii) lim(ii)], 'color', c, 'linewidth', 2);
end

if isempty(ols)==1 || val(1)==0
    h1 = [];
else
    h1 = plot(x, ols, '+', 'color', c, 'linewidth', 2);
end
h2 = plot([x-0.25 x-0.25], [lim(2) lim(4)],  'color', c, 'linewidth', 2);
h3 = plot([x+0.25 x+0.25], [lim(2) lim(4)], 'color', c, 'linewidth', 2);
h4 = plot([x x], [lim(4) lim(5)], 'color', c, 'linewidth', 2);
h5 = plot([x x], [lim(1) lim(2)], 'color', c, 'linewidth', 2);
h = [h' h1' h2' h3' h4' h5']';

plot([x-0.25 x+0.25], [lim1(2) lim1(2)], 'color', c, 'linewidth', 2);

end

