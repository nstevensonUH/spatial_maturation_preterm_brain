function fv1 = get_features(datm1, art1)



fs = 256;
A = size(datm1); %B = size(datm2);
olap = 900*fs;M = 60*60*fs;
block_no = floor(A(2)./olap)-3; if block_no<1; block_no=1; end
fv1 = cell(1, block_no); 
if block_no==1
    if length(datm1)> 60*60*fs; M = 60*60*fs; else; M = floor(length(datm1)/fs)*fs; end
    if sum(sum(art1(:,1:M))) < 0.5*A(1)*M
        fvz1 = extract_features_new_optim_v1(datm1(:, 1:M), fs, art1(:, 1:M));
        fv1{1} = fvz1;     
    else
        fv1{1} = NaN*ones(A(1),49);     
    end  
else
   for z1 = 1:block_no
       r1 = (z1-1)*olap+1; r2 = r1+3600*fs-1;
       if sum(sum(art1(:,r1:r2))) < 0.5*A(1)*M
           fvz1 = extract_features_new_optim_v1(datm1(:, r1:r2), fs, art1(:, r1:r2)); 
           fv1{z1} = fvz1;     
       else
           fv1{z1} = NaN*ones(A(1),49);     
       end
   end
end

