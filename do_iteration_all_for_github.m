function [r, ca, cb, outs1, outs2, err0, err1, err2, err3, err4, val_ref, ppx, out1, out2, out4, df] = do_iteration_all_for_github(fv2, pma2, id2, demos, MyAtlas, vlsx, fref)

zz=16;
out1 = cell(58); out2 = out1;  out4 = out1;
err1 = zeros(58); err2 = err1; err3 = err1; err4 = err1;
ppx = zeros(1,58);
for ii = 1:58
    ii
      Z1 = size(fv2); Z2 = size(fv2{1});
      fts1x = zeros(Z1(2), Z2(2));
        for z1 = 1:Z1(2)
            fts1x(z1,:) = fv2{z1}(ii,:);
        end
        pma1 = pma2; dem1 = demos;
        fts1x = fts1x(:, fref);
        ppx(ii) = corr(fts1x(:,9), pma1');
        % INITIAL FILTER STAGE
        D = size(fts1x); y = pma1'; rfs = zeros(D(2), D(1));
        for z1 = 1:D(2)
            x = fts1x(:, z1);
            B = robustfit(y, x);
            [~, dum2] = rmoutliers(x-(B(1)+B(2)*y));
            rfs(z1,:) = dum2';
        end
        nr = find(sum(rfs)<2);
        pma1x = pma1(nr);
        fts2x = fts1x(nr, :);
        id2x = dem1(nr,1); id2x(isnan(id2x)) = 2240;   

        fts2z = fts2x;
        pma1z = pma1x;
        id2z = id2x;
        
        [fts2x, pma1x, id2x] = augment_data_col(fts2x, pma1x', id2x, 9);

        responseScale = iqr(pma1);
        boxConstraint = responseScale/1.349;
        epsilon = responseScale/13.49;
        Mdl = fitrsvm(fts2x, pma1x, 'KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', zz, 'BoxConstraint', boxConstraint, 'Epsilon', epsilon, 'Standardize', true);  

        [X, Ydxc] = do_loso_model_col_v8(fts2x, pma1x, id2x, 1, zz, max(id2z));
        if vlsx==1
            Bxc = polyfit(X, Ydxc-X, 1);
            Ydx = Ydxc-polyval(Bxc, X);    
        end
                
    for kk = 1:58
        fts1 = zeros(Z1(2), Z2(2));
        for z1 = 1:length(fv2)
            fts1(z1,:) = fv2{z1}(kk,:);
        end
        pma1 = pma2; id3 = id2; rec3 = 1:length(fv2);
        fts1 = fts1(:,fref);
        pma1 = pma1(~isnan(sum(fts1'))); id3 = id3(~isnan(sum(fts1'))); rec3 = rec3(~isnan(sum(fts1')));
        fts1 = fts1(~isnan(sum(fts1')),:);
        D = size(fts1); y = pma1'; rfs = zeros(D(2), D(1));
        for z1 = 1:D(2)
            x = fts1(:, z1);
            B = robustfit(y, x);
            [~, dum2] = rmoutliers(x-(B(1)+B(2)*y));
            rfs(z1,:) = dum2';
        end
        nrx = find(sum(rfs)<2);
        pma1y = pma1(nrx);
        fts2y = fts1(nrx, :);
        %demy = dem3(nrx,:);             
        Ydc = predict(Mdl, fts2y);
        cad = pma1y';
        
        [~, ia, ib] = intersect(nr, nrx);
        
        if vlsx == 1
            
        Yd = Ydc-polyval(Bxc, cad);
        out1{ii,kk} = Ydx(ia); out2{ii,kk} = X(ia); 
        out4{ii,kk} = Yd(ib)';
        err0(ii,kk) = wbias(Ydx(ia), Yd(ib)');
        err1(ii,kk) = mean((Ydx(ia) - Yd(ib)'));
        %err0(ii,kk) = wbias(Yd', cad');
        %err1(ii,kk) = mean((Yd-cad));
        err2(ii,kk) = std(Yd-cad);
        err3(ii,kk) = corr(Yd, cad);
        err4(ii,kk) = mean(abs(Yd-cad));
        
        else
            
        out1{ii,kk} = Ydxc(ia); out2{ii,kk} = X(ia); 
        out4{ii,kk} = Ydc(ib)';
        err0(ii,kk) = wbias(Ydxc(ia), Ydc(ib)');
        err1(ii,kk) = mean((Ydxc(ia) - Ydc(ib)'));
        %err0(ii,kk) = wbias(Ydc', cad');
        %err1(ii,kk) = mean((Ydc-cad));
        err2(ii,kk) = std(Ydc-cad);
        err3(ii,kk) = corr(Ydc, cad);
        err4(ii,kk) = mean(abs(Ydc-cad));
            
        end

    end    

    if vlsx == 1
    
    err0(ii,ii) = wbias(Ydx, X);
    err1(ii,ii) = 0;
    err2(ii,ii) = std(Ydx-X);
    err3(ii,ii) = corr(Ydx',X');
    err4(ii,ii) = mean(abs(Ydx-X));
    out1{ii,ii} = Ydx; out2{ii,ii} = X; 
    out4{ii,ii} = nr;
    
    else
        
    err0(ii,ii) = wbias(Ydxc, X);
    err1(ii,ii) = 0;
    err2(ii,ii) = std(Ydxc-X);
    err3(ii,ii) = corr(Ydxc',X');
    err4(ii,ii) = mean(abs(Ydxc-X));
    out1{ii,ii} = Ydxc; out2{ii,ii} = X; 
    out4{ii,ii} = nr;
    
    end


end

dum = MyAtlas.Centroids;
idx = zeros(1,58);
idx([19 20 22 23 25 28 29 48 49 51 52 54 57 58]) = 2;  %sources at depth
val_ref = find(idx==0);
dum(:,2) = abs(dum(:,2));
dum  = dum(val_ref,:);
sar = MyAtlas.SArank(val_ref);

outs1 = cell(4,2); outs2 = outs1;
pad = zeros(1, length(val_ref)); par = pad; pas = par; 
padw = pad;
for ii = 1:length(val_ref)
    rf = 1:length(val_ref);
    rf(ii)=0; rf1 = find(rf>0);
    
    padw(ii) = median(err0(val_ref(ii), val_ref(rf1)));
    pad(ii) = median(err1(val_ref(ii), val_ref(rf1)));
    par(ii) = median(err3(val_ref(ii), val_ref(rf1)));    
    pas(ii) = median(err2(val_ref(ii), val_ref(rf1)));    
end
outs1{1,1} = pad; outs1{1,2} = dum;
[r1d, ~, c1d, c2d] = corrcoef(dum(:,1), padw);
r0(1) = r1d(1,2); c0a(1) = c1d(1,2); c0b(1) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(dum(:,2), padw);
r0(2) = r1d(1,2); c0a(2) = c1d(1,2); c0b(2) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(dum(:,3), padw);
r0(3) = r1d(1,2); c0a(3) = c1d(1,2); c0b(3) = c2d(1,2);


[r1d, ~, c1d, c2d] = corrcoef(dum(:,1), pad);
r1(1) = r1d(1,2); c1a(1) = c1d(1,2); c1b(1) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(dum(:,2), pad);
r1(2) = r1d(1,2); c1a(2) = c1d(1,2); c1b(2) = c2d(1,2);

[r1d, ~, c1d, c2d] = corrcoef(dum(:,3), pad);
r1(3) = r1d(1,2); c1a(3) = c1d(1,2); c1b(3) = c2d(1,2);

rng(7)
rf = zeros(1000,length(val_ref));
for ii = 1:1000
    rf(ii,:) = randsample(length(val_ref), length(val_ref), true);
end
df1 = zeros(1000,4);
for ii = 1:1000
    B = polyfit(dum(rf(ii,:),1), pad(rf(ii,:)), 1);
    df1(ii,1) = polyval(B,max(dum(:,1))) - polyval(B, min(dum(:,1)));   
    B = polyfit(dum(rf(ii,:),2), pad(rf(ii,:)), 1);
    df1(ii,2) = polyval(B,max(dum(:,2))) - polyval(B, min(dum(:,2)));   
    B = polyfit(dum(rf(ii,:),3), pad(rf(ii,:)), 1);
    df1(ii,3) = polyval(B,max(dum(:,3))) - polyval(B, min(dum(:,3)));   
    B = polyfit(sar(rf(ii,:)), pad(rf(ii,:)), 1);
    df1(ii,4) = polyval(B,max(sar)) - polyval(B, min(sar));   
end

df = [mean(df1)' quantile(df1, [0.025 0.975])'];

outs1{2,1} = pad; outs1{2,2} = sar;
[r1d, ~, c1d, c2d] = corrcoef(sar, padw);
r0(4) = r1d(1,2); c0a(4) = c1d(1,2); c0b(4) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(sar, pad);
r1(4) = r1d(1,2); c1a(4) = c1d(1,2); c1b(4) = c2d(1,2);


outs2{1,1} = par; outs2{1,2} = dum;
outs2{2,1} = par; outs2{2,2} = sar;
[r1d, ~, c1d, c2d] = corrcoef(dum(:,1), par);
r2(1) = r1d(1,2); c2a(1) = c1d(1,2); c2b(1) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(dum(:,2), par);
r2(2) = r1d(1,2); c2a(2) = c1d(1,2); c2b(2) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(dum(:,3), par);
r2(3) = r1d(1,2); c2a(3) = c1d(1,2); c2b(3) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(sar, par);
r2(4) = r1d(1,2); c2a(4) = c1d(1,2); c2b(4) = c2d(1,2);

[r1d, ~, c1d, c2d] = corrcoef(dum(:,1), pas);
r3(1) = r1d(1,2); c3a(1) = c1d(1,2); c3b(1) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(dum(:,2), pas);
r3(2) = r1d(1,2); c3a(2) = c1d(1,2); c3b(2) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(dum(:,3), pas);
r3(3) = r1d(1,2); c3a(3) = c1d(1,2); c3b(3) = c2d(1,2);
[r1d, ~, c1d, c2d] = corrcoef(sar, pas);
r3(4) = r1d(1,2); c3a(4) = c1d(1,2); c3b(4) = c2d(1,2);

r = [r0 ; r1 ; r2 ; r3];
ca = [c0a ; c1a ; c2a ; c3a];
cb = [c0b ; c1b ; c2b ; c3b];

