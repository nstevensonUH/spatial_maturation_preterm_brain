function [new_feature_set, ages, ids] = augment_data_col(f1, pm1, idx, MM)
%
%  For features with a correlation with age
%  Needs adjustment if covariance is age dependent, which can probably only
%  be assessed if the data is big
%
%

NN = size(f1);
stp = (max(pm1)-min(pm1))./MM;
rnge = [min(pm1):stp:max(pm1)-stp ; min(pm1)+stp:stp:max(pm1)];
rnge(2,MM) = max(pm1)+1;
hs = zeros(1,MM);
for ii = 1:MM
    hs(ii) = length(find(pm1>=rnge(1,ii) & pm1<rnge(2,ii)+eps));
end
valx = max(hs)-hs;

% use covariance you dingbat so generate mvrndn for each sample
dc = f1;
pm = pm1; M = length(pm);
PP = cell(1,NN(2));
dum = zeros(NN(2),M); idy = dum;
for kk = 1:NN(2)
    aic = zeros(1,3);
    ndc = dc(:,kk);
    val = zeros(1, MM);
    for ii = 1:MM
        rf = find(pm1>=rnge(1,ii) & pm1<rnge(2,ii)+eps);
        val(ii) = median(rmoutliers(ndc(rf)), 'omitnan');
    end
    for ii = 1:4
        P = polyfit(mean(rnge), val, ii-1);
        Y = polyval(P, pm);
        nu = (M-ii-1);
        aic(ii) = (sum(((dc(:,kk)-Y).^2)))./nu;
    end    
    val = find(aic==min(aic));
    PP{kk} = polyfit(pm, dc(:,kk), val-1);
    Y = polyval(PP{kk}, pm);
    dum(kk,:) = (ndc-Y)';
    [~, idy(kk,:)] = rmoutliers(dum(kk,:));
end

sig = cov(dum(:,find(sum(idy)==0))');
rng(3)
% Augment data with synthetic values
nagx = []; nftx = []; idy = []; c1 = (max(idx)/10)+1;
for ii = 1:MM
    for jj = 1:valx(ii)
     agx = rand(1)*stp+rnge(1,ii);
     mx = zeros(1,NN(2));
     for kk = 1:NN(2)
        mx(kk) = polyval(PP{kk}, agx);
     end
     out = mvnrnd(mx, sig, 1);
     nagx = [nagx ; agx];
     nftx = [nftx ; out];
     idy = [idy ; c1*10]; c1 = c1+1;
    end
end

new_feature_set = [f1 ; nftx];
ages = [pm1 ; nagx];
ids = [idx ; idy];
