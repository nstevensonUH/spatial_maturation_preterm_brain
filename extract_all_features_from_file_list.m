function fvx1= extract_all_features_from_file_list(aa)
% This function performs source space reconstruction and then extracts a 
% list of features from a series of pre-process EDF files defined
% by variable aa and saved in mat files at a sampling frequency of 256Hz
%
% Inputs: aa -> a cell array containing a list of files, these input files
% are EEG with artefact removed or identified of approximately 1h in length
% with 125 channels
%
% Outputs: fvx1 -> a cell array containing 49 features calculated on each
% channel/region/parcel per recording/filename, i.e. 58 x 49 matrix per
% cell entry
% 
%
%
% Dependeincies - requires mat files CollapseOperator.mat,
%                                    InverseOperator.mat, and 
%                                    MyAtlas_n58.mat (on Github in supporting_mat_files folder)
%
% dependencies on Github in supporting_m_files folder
%  - get_features.m
%       - extract_features_new_optim_v1.m
%           - detector_per_channel_palmu_adj.m
%               - nlin_energy.m
%               - precess_ba_1ch.m
%           - estimate_mse_och_fast.m
%               - SampEn.m
%
%
% Nathan Stevenson 
% QIMR Berghofer

% These files are on the github repository, add extra directory information
% at the front of filenames if not processing in the same directory as this
% files or the directory containing thises files is not added to the path
load('CollapseOperator.mat')
load('InverseOperator.mat')
load('MyAtlas_n58.mat')


N = 1024; fs1 = 256;
CollapseOperator = repmat(CollapseOperator, 1, N);
fvx1 = cell(1,length(aa));
for zz = 1:length(aa)

    filename = aa{zz};
    load(filename)
    dumx = dumy;        
    clear dumy

    % DO SOURCE SPACE RECONSTRUCTION

    dat5 = zeros(58,length(dumx));
    block_no = floor(length(dumx)/N); Np = length(MyAtlas.Parcels); % number of parcels
    for ii = 1:block_no
        r1 = (ii-1)*N+1; r2 = r1+N-1;
        dumy = dumx(:, r1:r2);
        % source signals  
        src_buf = (InverseOperator * dumy).*CollapseOperator;     % [sources x samples]
        % parcel/cortical signals  
        parcel_buf = zeros(Np, N);
        for jj = 1:Np
           parcel_buf(jj, :) = mean(src_buf(MyAtlas.Parcels{jj, 1}, :));     % [parcels x samples] 
        end
        dat5(:,r1:r2) = parcel_buf;
    end
    clear src_buf parcel_buf

    % LOOK FOR HIGH AMPLITUDE ARTEFACTS IN DATA - to ignore in future
    % processing

    block_no = floor(length(dat5)/fs1); artdum = zeros(58,block_no);
    for ii = 1:block_no
        r1 = (ii-1)*fs1+1; r2 = r1+fs1-1;
        for jj = 1:58
            dum = abs(hilbert(dat5(jj, r1:r2)));
            artdum(jj,ii) = max(max(dum));
        end
    end
    duma = rmoutliers(artdum(:));
    amp_norm_X = median(duma);
    dat5 = dat5./amp_norm_X.*median(rmoutliers(abs(dumx(1,:))));
    clear dumx
    art5 = zeros(58,block_no);
    for ii = 1:block_no
        r1 = (ii-1)*fs1+1; r2 = r1+fs1-1;
        for jj = 1:58
            dum = abs(hilbert(dat5(jj, r1:r2)));
            art5(jj,ii) = max(max(dum));
        end
    end
    art5(art5<200)=0; art5(art5>=200)=1; % changed to 200 due to filtering first
    for ii = 1:58
       dum = art5(ii,:);
       r1 = find(diff([0 dum 0])==1);
       r2 = find(diff([0 dum 0])==-1);
       r1 = r1-2; r2 = r2+2;
       r1(r1<1)=1;
       r2(r2>length(dum))=length(dum);
       dum1 = zeros(size(dum));
       for jj = 1:length(r1)
           dum1(r1(jj):r2(jj))=1;
       end
       art5(ii,:)=dum1;
    end
    
    % EXTRACT FEATURES FROM SSR DATA
    fvx1{zz} = get_features(dat5, art5);
    
    clear dat5 art5 

end