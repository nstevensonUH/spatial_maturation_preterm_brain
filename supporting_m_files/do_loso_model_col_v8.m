function [X, Y] = do_loso_model_col_v8(fts2, ca2, id2, flag, zz, max_v)


M = 10;
yr = divide_data_col(id2', ca2', M);
rng(13)
% 5 fold CV here
Y = zeros(1, find(id2==max_v, 1, 'last')); X = Y;
for ii = 1:M
    
    ref1 = []; for jj = 1:length(yr{ii}); ref1 = [ref1 ; find(id2==yr{ii}(jj))];  end
    didx = 1:length(id2);
    didx(ref1) = 0;
    ref2 = find(didx~=0);
     
    pma1 = ca2(ref2);

    switch flag
        case 1
        responseScale = iqr(pma1);
        boxConstraint = responseScale/1.349;
        epsilon = responseScale/13.49;
        Mdl = fitrsvm(fts2(ref2,:), pma1', 'KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', zz, 'BoxConstraint', boxConstraint, 'Epsilon', epsilon, 'Standardize', true);  
        case 3
        responseScale = iqr(pma1);
        boxConstraint = responseScale/1.349;
        epsilon = responseScale/13.49;
        Mdl = fitrsvm(fts2(ref2,:), pma1', 'KernelFunction', 'linear', 'BoxConstraint', boxConstraint, 'Epsilon', epsilon, 'Standardize', true); 
        otherwise
        Mdl = fitrgp(fts2(ref2,:), pma1', 'BasisFunction', 'constant', 'KernelFunction', 'rationalquadratic', 'Standardize', true);
    end
    
    ref3 = find(id2(ref1)<=max_v);
    dum = predict(Mdl, fts2(ref1(ref3),:));
    Y(ref1(ref3)) = dum;
    X(ref1(ref3)) = ca2(ref1(ref3));
end


% pma1 = ca2;
% nrx = find(pma1>38);
% rngs = [32 36 ; 36 38]; nn = 50;
% for pp = 1:2
%      rxx = find(pma1>rngs(1,pp) & pma1<=rngs(2,pp));
%      if nn<length(rxx); yrx = randsample(length(rxx), nn1); else;  yrx = 1:length(rxx); end
%      nrx = [nrx rxx(yrx)];
% end
% 
% res = Y(nrx) - pma1(nrx);
% Bxc = regress(res', [pma1(nrx)' ones(length(nrx),1)]);
% Ydxc = Y-(Bxc(1)*pma1+Bxc(2));
        

