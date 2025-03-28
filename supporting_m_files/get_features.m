function fv1 = get_features(datm1, art1)



fs = 256;
A = size(datm1); B = size(datm2);
olap = 900*fs;
block_no = floor(A(2)./olap)-3; if block_no<1; block_no=1; end
ar1 = sum(art1); ar1(ar1<20)=0; ar1(ar1>=20)=1; 
fv1 = cell(1, block_no); 
if block_no==1
    if length(datm1)> 60*60*fs; M = 60*60*fs; else; M = floor(length(datm1)/fs)*fs; end
    if sum(ar1(1:M/fs)) < 0.5*(length(datm1(:,1:M))/fs)
        fvz1 = extract_features_new_optim_v1(datm1(:, 1:M), fs, art1(:, 1:floor(M/fs))); % burst_measures_new_plus_ema_v1_columbia(datm1(:, 1:M), fs, art1(:,1:floor(M/fs)));
        fv1{1} = fvz1;     
    else
        fv1{1} = NaN*ones(A(1),49);     
    end  
else
   for z1 = 1:block_no
       r1 = (z1-1)*olap+1; r2 = r1+3600*fs-1;
       if sum(ar1(ceil(r1/fs):r2/fs)) < 0.5*length(r1:r2)/fs
           fvz1 = extract_features_new_optim_v1(datm1(:, 1:M), fs, art1(:, 1:floor(M/fs))); % burst_measures_new_plus_ema_v1_columbia(datm1(:, r1:r2), fs, art1(:,ceil(r1/fs):r2/fs));
           fv1{z1} = fvz1;     
       else
           fv1{z1} = NaN*ones(A(1),49);     
       end
   end
end

