function out = extract_features_new_optim_v1(data, fs1, art)
%addpath(genpath('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\BRM3'))
%addpath('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\burst_analysis\for_Sampsa\Innsbruck\')
%addpath('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\burst_analysis\code\final')
%fname = 'L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Helsinki\ElkeG\ANONYM_20210417_RAW.brm'

% EXTRACT FEATURES (1h epochs, 15 minute overlap)
fs2 = 64;
dat64 = resample(data', fs2, fs1)';
a64 = resample(art', fs2, fs1)'; a64(a64>0.5)=1; a64(a64<1)=0;
durbins=[1 8 ; 8 16 ; 16 32 ; 32 64 ; 64 128 ; 128 256 ; 1 8192]*fs1/fs2;
ob = [0.5 2 ; 2 4 ; 4 8 ; 8 13 ; 13 32];
dlim = 5; % CHANGE TO xlim burst fit, potential limit based on CDF analysis (35ms at the moment - half a fast spike)
A = size(data);
feat = zeros(A(1),49); valid = zeros(A(1),1);
    for ch = 1:A(1)
        try
         dat = dat64(ch, :);
         dat1 = data(ch, :);
         adat = a64(ch, :);
         adat1 = art(ch, :);
         if sum(adat)<0.5*length(adat) & sum(abs(dat1))>0
         % AMPLITUDE
         dum = quantile(abs(hilbert(dat)), [0.05 0.25 0.5 0.75 0.95]);
         feat(ch,1) = dum(3);
         feat(ch, 2:5) = dum([1 2 4 5])';
         % BURST DURATION & INTERVAL
         ba = detector_per_channel_palmu_adj(dat, fs2, adat); % Using Kirsi's burst detector. % CHECK FOR SPEED
         r1z = find(diff([0 ba 0]) == 1);
         r2z = find(diff([0 ba 0]) == -1);
         ibis = r1z(2:end)-r1z(1:end-1);
         bdurs = r2z-r1z;

         a1z = find(diff([0 adat 0]) == 1);
         a2z = find(diff([0 adat 0]) == -1);
         refx = 1:length(bdurs);
         for pp = 1:length(a1z)
             ref1 = [find(abs(r1z-a1z(pp)) == min(abs(r1z-a1z(pp))),1) find(abs(r2z-a1z(pp)) == min(abs(r2z-a1z(pp))),1)];
             ref2 = [find(abs(r1z-a2z(pp)) == min(abs(r1z-a2z(pp))),1) find(abs(r2z-a2z(pp)) == min(abs(r2z-a2z(pp))),1)];
             val = [min([ref1 ref2])-1 max([ref1 ref2])+1];    
             if val(1)<1; val(1)=1; end
             if val(end)>length(bdurs); val(end)=length(bdurs); end
             refx(val(1):val(2))=0;
         end
         feat(ch, 6:8) = quantile(ibis(refx(1:end-1)>0), [0.05 0.5 0.95])/fs2;
         feat(ch, 9:11) = quantile(bdurs(refx>0), [0.05 0.5 0.95])/fs2;

         % SKW / KURT BURSTS 
        amp = abs(hilbert(dat1)).^2;
        amp(adat==1)=0;
        amp = conv(amp, [1 1 1 1 1]/5, 'same');
        th1 = quantile(amp, 100); r1q = zeros(1,length(th1));
        for zz = 1:length(th1)
            dum = zeros(1, length(amp));
            dum(amp<th1(zz)) = 0;
            dum(amp>=th1(zz)) = 1;      
            r1x = find(diff([0 dum 0]) == 1);
            r2x = find(diff([0 dum 0]) == -1);
            tst2 = r2x-r1x;
            rf = find(tst2<=dlim);
            for z2 = 1:length(rf)
                dum(r1x(rf(z2)):r2x(rf(z2)))=0;
            end
            r1q(zz) = length(find(diff([0 dum 0])==1));
        end
        th = th1(find(r1q==max(r1q),1)); 
        dum(amp<th) = 0;
        dum(amp>=th) = 1;
        r1x = find(diff([0 dum 0]) == 1);
        r2x = find(diff([0 dum 0]) == -1);
        tst2 = r2x-r1x;
        r1x = r1x(tst2>dlim);
        r2x = r2x(tst2>dlim);
        tst2 = tst2(tst2>dlim);
        M = ceil(2*mean(tst2)*10);
        sk = zeros(1,length(durbins)); kt = sk;
        if isfinite(M)==1
        for qq = 1:length(durbins)
            bav = zeros(1,M);
            r1y = r1x(tst2 > durbins(qq,1) & tst2 <= durbins(qq,2));
            r2y = r2x(tst2 > durbins(qq,1) & tst2 <= durbins(qq,2));
            for kk = 1:length(r1y)
               dum =  amp(r1y(kk):r2y(kk)-1)-th; 
               xx = linspace(1,length(dum),M);
               pp = pchip(1:length(dum), dum, xx); 
               pp = pp./sum(pp);
               bav = bav+pp; % Average shape
            end
            % 4 moments of burst shape only need skewness and kurtosis
            bav1 = (bav-min(bav)); bav1 = bav1./sum(bav1);
            xx = linspace(0,1,M);
            mn = sum(xx.*bav1);
            sd = sqrt(sum((xx-mn).^2.*bav1));
            sk(qq) = sum((xx-mn).^3.*bav1)./sd^3;
            kt(qq) = sum((xx-mn).^4.*bav1)./sd^4-3;
        end
        end
        feat(ch, 12:18) = sk;
        feat(ch, 19:25) = kt;

        bd1 = zeros(1,length(r1x)); ba1 = bd1;
        for kk= 1:length(r1x)
            bd1(kk) = length(r1x(kk):r2x(kk)-1)/fs2;
            ba1(kk) = trapz(amp(r1x(kk):r2x(kk)-1)-th)/fs2; % Area Trapezoidal
        end

        B1 = polyfit(log(bd1), log(ba1), 1);
        B2 = polyfit(durbins(:,1)'./fs2, sk, 1);
        B3 = polyfit(durbins(:,1)'./fs2, kt, 1); 
        arx = resample(adat1, fs1, 1); arx(arx>0.5)=1; arx(arx<=0.5)=0;
        SC = calculateSC_NS_pch(dat1, fs1, arx);
        bdm = mean(tst2)/fs2;
        bdst = std(tst2)/fs2;

        feat(ch, 26:34) = [B1 B2 B3 SC bdm bdst];

        % power laws 
        dumb = zeros(1,3);
        MM = ceil(sqrt(length(bd1))); cv = floor(MM/12);
        [b, n] = ucdf_log(bd1, MM);
        y1 = log10(b)'; x1 = n(isfinite(y1))'; y1 = y1(isfinite(y1));
        res = ones(length(y1), length(y1));%mean((y1-mean(y1)).^2)*ones(length(y1), length(y1));
        for z1 = cv+1:length(y1)
            if z1+cv+16 > length(y1)-cv; len = length(y1)-cv; else; len = z1+cv+16; end
            for z2 = z1+cv:len
                B = regress(y1(z1:z2), [x1(z1:z2) ones(size(x1(z1:z2)))]);
                z11 = z1-cv; z22 = z2+cv; ids = [z11:z1-1 z2+1:z22];
                res1 = y1(ids)' - (B(1)*x1(ids)'+B(2)*ones(1,length(x1(ids))));
                res(z1,z2) = mean(res1.^2)/length(z1:z2)^2;            
            end
        end
        [q1, q2] = find(res==min(min(res)));
        z11 = q1-cv; z22 = q2+cv; ids = z11:z22;
        B = regress(y1(ids), [x1(ids) ones(size(x1(ids)))]);            
        dumb(1) = B(1);
        y1 = log10(b)'; x1 = n(isfinite(y1))'; y1 = y1(isfinite(y1));
        y1(ids)=NaN; x1(ids)=NaN; y1 = y1(isnan(y1)==0); x1 = x1(isnan(x1)==0);
        y2 = y1; x2 = x1;
        res = ones(length(y1), length(y1));%mean((y1-mean(y1)).^2)*ones(length(y1), length(y1));
        for z1 = cv+1:length(y1)
            if z1+cv+16 > length(y1)-cv; len = length(y1)-cv; else; len = z1+cv+16; end
            for z2 = z1+cv:len
                B = regress(y1(z1:z2), [x1(z1:z2) ones(size(x1(z1:z2)))]);
                z11 = z1-cv; z22 = z2+cv; ids = [z11:z1-1 z2+1:z22];
                res1 = y1(ids)' - (B(1)*x1(ids)'+B(2)*ones(1,length(x1(ids))));
                res(z1,z2) = mean(res1.^2)/length(z1:z2)^2;            
            end
        end
        [q1, q2] = find(res==min(min(res)),1);
        z11 = q1-cv; z22 = q2+cv; ids = z11:z22;
        B = regress(y1(ids), [x1(ids) ones(size(x1(ids)))]); 
        dumb(2) = B(1);
        try
        y1 = y2; x1= x2;
        y1(ids)=NaN; x1(ids)=NaN; y1 = y1(isnan(y1)==0); x1 = x1(isnan(x1)==0);
        res = ones(length(y1), length(y1));%mean((y1-mean(y1)).^2)*ones(length(y1), length(y1));
        for z1 = cv+1:length(y1)
            if z1+cv+16 > length(y1)-cv; len = length(y1)-cv; else; len = z1+cv+16; end
            for z2 = z1+cv:len
                B = regress(y1(z1:z2), [x1(z1:z2) ones(size(x1(z1:z2)))]);
                z11 = z1-cv; z22 = z2+cv; ids = [z11:z1-1 z2+1:z22];
                res1 = y1(ids)' - (B(1)*x1(ids)'+B(2)*ones(1,length(x1(ids))));
                res(z1,z2) = mean(res1.^2)/length(z1:z2)^2;            
            end
        end
        [q1, q2] = find(res==min(min(res)),1);
        z11 = q1-cv; z22 = q2+cv; ids = z11:z22;
        B = regress(y1(ids), [x1(ids) ones(size(x1(ids)))]); 
        dumb(3) = B(1);
        feat(ch, 35:37) = sort(dumb, 'descend');
        catch
        feat(ch, 35:37) = sort(dumb, 'descend'); 
        end

        clear dumb x1 y1 x2 y2
        dumb = zeros(1,3);
        MM = ceil(sqrt(length(ba1))); cv = floor(MM/12);
        [b, n] = ucdf_log(ba1, MM);
        y1 = log10(b)'; x1 = n(isfinite(y1))'; y1 = y1(isfinite(y1));
        res = ones(length(y1), length(y1));%mean((y1-mean(y1)).^2)*ones(length(y1), length(y1));
        for z1 = cv+1:length(y1)
            if z1+cv+16 > length(y1)-cv; len = length(y1)-cv; else; len = z1+cv+16; end
            for z2 = z1+cv:len
                B = regress(y1(z1:z2), [x1(z1:z2) ones(size(x1(z1:z2)))]);
                z11 = z1-cv; z22 = z2+cv; ids = [z11:z1-1 z2+1:z22];
                res1 = y1(ids)' - (B(1)*x1(ids)'+B(2)*ones(1,length(x1(ids))));
                res(z1,z2) = mean(res1.^2)/length(z1:z2);            
            end
        end
        [q1, q2] = find(res==min(min(res)));
        z11 = q1-cv; z22 = q2+cv; ids = z11:z22;
        B = regress(y1(ids), [x1(ids) ones(size(x1(ids)))]);            
        dumb(1) = B(1);
        y1 = log10(b)'; x1 = n(isfinite(y1))'; y1 = y1(isfinite(y1));
        y1(ids)=NaN; x1(ids)=NaN; y1 = y1(isnan(y1)==0); x1 = x1(isnan(x1)==0);
        y2 = y1; x2 = x1;
        res = ones(length(y1), length(y1));%mean((y1-mean(y1)).^2)*ones(length(y1), length(y1));
        for z1 = cv+1:length(y1)
            if z1+cv+16 > length(y1)-cv; len = length(y1)-cv; else; len = z1+cv+16; end
            for z2 = z1+cv:len
                B = regress(y1(z1:z2), [x1(z1:z2) ones(size(x1(z1:z2)))]);
                z11 = z1-cv; z22 = z2+cv; ids = [z11:z1-1 z2+1:z22];
                res1 = y1(ids)' - (B(1)*x1(ids)'+B(2)*ones(1,length(x1(ids))));
                res(z1,z2) = mean(res1.^2)/length(z1:z2);            
            end
        end
        [q1, q2] = find(res==min(min(res)),1);
        z11 = q1-cv; z22 = q2+cv; ids = z11:z22;
        B = regress(y1(ids), [x1(ids) ones(size(x1(ids)))]); 
        dumb(2) = B(1);
        try
        y1 = y2; x1= x2;
        y1(ids)=NaN; x1(ids)=NaN; y1 = y1(isnan(y1)==0); x1 = x1(isnan(x1)==0);
        res = ones(length(y1), length(y1));%mean((y1-mean(y1)).^2)*ones(length(y1), length(y1));
        for z1 = cv+1:length(y1)
            if z1+cv+16 > length(y1)-cv; len = length(y1)-cv; else; len = z1+cv+16; end
            for z2 = z1+cv:len
                B = regress(y1(z1:z2), [x1(z1:z2) ones(size(x1(z1:z2)))]);
                z11 = z1-cv; z22 = z2+cv; ids = [z11:z1-1 z2+1:z22];
                res1 = y1(ids)' - (B(1)*x1(ids)'+B(2)*ones(1,length(x1(ids))));
                res(z1,z2) = mean(res1.^2)/length(z1:z2);            
            end
        end
        [q1, q2] = find(res==min(min(res)),1);
        z11 = q1-cv; z22 = q2+cv; ids = z11:z22;
        B = regress(y1(ids), [x1(ids) ones(size(x1(ids)))]); 
        dumb(3) = B(1);
        feat(ch, 38:40) = sort(dumb, 'descend');
        catch
        feat(ch, 38:40) = sort(dumb, 'descend'); 
        end

        % Sample Entropy
        feat(ch,41) = estimate_mse_pch_fast(dat1, fs1, arx);

        % SPECTRAL ANALYSIS, POWER/SLOPE/SPECTRAL EDGE (90%)
        N1 = length(dat1);
        f = 0:fs1/N1:fs1-1/N1;
        BB = size(ob);
        % Estimate PSD - periodogram across all channels
        dta = dat1;
        dta(arx==1) = 0;
        dta = dta-mean(dta);
        N = sum(arx==0);
        DTA = 1/N*abs(fft(dta)).^2; % periodogram (slighlty modified due to artefact)
        DTA = 2*DTA(:,1:ceil(N1/2)); % still works as DTA(1) will be zero
        fref1 = find(f>=min(min(ob)) & f<max(max(ob))); % full band of interest
        tp = sum(DTA(fref1));
        feat(ch,42) = log10(tp);
        % Estimate relative spectral powers in each band
        ff = zeros(1,BB(1));
        for z2 = 1:BB(1)
           fref1 = find(f>=ob(z2,1) & f<ob(z2, 2)); 
           ff(z2) = sum(DTA(fref1))./tp;
        end
        feat(ch, 43:47) = ff*100;

        [Pxx, f] = pwelch(dta, hamming(2^16), 2^15, 2^16, fs1);
        fr = 2.^[log2(min(min(ob))+1.5):0.1:log2(16)]; P = zeros(1,length(fr)-1);
        for cc = 1:length(fr)-1
           rx = find(f>=fr(cc) & f<(fr(cc+1)));
           P(cc) = mean(Pxx(rx)); 
        end  
       outliers = 1; logf = log2(fr); logP = log2(P); idx = zeros(1,length(P));
       while outliers ~= 0
            x1 = logf(idx==0); y1 = logP(idx==0);
            B = regress(y1', [x1' ones(size(x1'))]);
            res = y1 - (B(1)*x1+B(2)*ones(1,length(x1)));
            [~, idx] = rmoutliers(res);
            outliers = sum(idx);
       end
       feat(ch,48) = B(1);
       fdum = f(f>=0.5 & f<=32);
       frx = find(cumsum(Pxx(f>=0.5 & f<=32))./sum(Pxx(f>=0.5 & f<=32))>0.9, 1); 
       feat(ch,49) = fdum(frx); 
       else
             valid(ch)=1;
         end
       catch
           valid(ch) = 1;
       end
    end

out = cell(1,2);
out{1} = feat;
out{2} = valid;

end    

% flist{1} = 'Median Amplitude Envelope'; %Amplitude Envelope is estimate using the analytic associate of a signal (via Hilbert transform)
% flist{2} = '5th Percentile Amplitude Envelope';
% flist{3} = '25th Percentile Amplitude Envelope';
% flist{4} = '75th Percentile Amplitude Envelope';
% flist{5} = '95th Percentile Amplitude Envelope';
% flist{6} = '5th Percentile Inter-burst Interval';
% flist{7} = '50th Percentile Inter-burst Interval';
% flist{8} = '95th Percentile Inter-burst Interval';
% flist{9} = '5th Percentile Burst Duration';
% flist{10} = '50th Percentile Burst Duration';
% flist{11} = '95th Percentile Burst Duration';
% flist{12} = 'Burst Skewness/Symmetry (15.625-125ms)';
% flist{13} = 'Burst Skewness/Symmetry (125-250ms)';
% flist{14} = 'Burst Skewness/Symmetry (250-500ms)';
% flist{15} = 'Burst Skewness/Symmetry (0.5-1s)';
% flist{16} = 'Burst Skewness/Symmetry (1-2s)';
% flist{17} = 'Burst Skewness/Symmetry (2-4s)';
% flist{18} = 'Burst Skewness/Symmetry (15.625ms-128s)';
% flist{19} = 'Burst Kurtosis/Sharpness (15.625-125ms)';
% flist{20} = 'Burst Kurtosis/Sharpness (125-250ms)';
% flist{21} = 'Burst Kurtosis/Sharpness (250-500ms)';
% flist{22} = 'Burst Kurtosis/Sharpness (0.5-1s)';
% flist{23} = 'Burst Kurtosis/Sharpness (1-2s)';
% flist{24} = 'Burst Kurtosis/Sharpness (2-4s)';
% flist{25} = 'Burst Kurtosis/Sharpness (15.625ms-128s)';
% flist{26} = 'Slope Burst Duration vs Burst Area';
% flist{27} = 'Intercept Burst Duration vs Burst Area';
% flist{28} = 'Slope Burst Duration vs Skewness';
% flist{29} = 'Intercept Burst Duration vs Skewness';
% flist{30} = 'Slope Burst Duration vs Kurtosis';
% flist{31} = 'Intercept Burst Duration vs Kurtosis';
% flist{32} = 'Suppression Curve';
% flist{33} = 'Mean Burst Duration';
% flist{34} = 'Standard Deviation Burst Durations';
% flist{35} = 'Power Law of Burst Duration Distribution alpha 1';
% flist{36} = 'Power Law of Burst Duration Distribution alpha 2';
% flist{37} = 'Power Law of Burst Duration Distribution alpha 3';
% flist{38} = 'Power Law of Burst Area Distribution alpha 1';
% flist{39} = 'Power Law of Burst Area Distribution alpha 2';
% flist{40} = 'Power Law of Burst Area Distribution alpha 3';
% flist{41} = 'Sample Entropy';
% flist{42} = 'Log PSD Energy (wideband)';
% flist{43} = 'Relative Delta 1 power';
% flist{44} = 'Relative Delta 2 power';
% flist{45} = 'Relative Theta power';
% flist{46} = 'Relative Alpha power';
% flist{47} = 'Relative Beta power';
% flist{48} = 'Spectral Slope';
% flist{49} = 'Spectral Edge Frequency (90th Percentile)';
% 
    