function bm = detector_per_channel_palmu_adj(eeg_data, Fs, artx)
% Kirsi's EEG 'burst' detector. Code based on the paper - Palmu et al,
% Physiol Meas 31:85-93, 2010.
%
% [t, bm] = detector_per_channel_palmu(eeg_data, Fs)
%
% Inputs:- a single channel of EEG (vector 1xN)
%              - sampling frequency of EEG
%
% Outputs - time vector (output is sampled at 16Hz - vector 1xM)
%                 - binary burst mask (1 - burst, 0 - not burst - vector 1xM)
%
% Nathan Stevenson
% August 2016
% University of Helsinki

% Initialise Parameters (feel free to tinker)
win = 1.5;
wc = 20;

load kirsi_BD_filters_256 % HP 6th order Elliptical fc = 10Hz. LP 1st order Butterworth fc = 0.5Hz

% Initialise Parameters
fs1 = 256; fs2 = 16;
dat = resample(eeg_data, fs1, Fs); % NB - must be at 256 as NLEO is different when when fs changes so parameters above will need to change is fs is changed.
arx = resample(artx, fs2, 1);
arx(arx<0.5)=0; arx(arx>=0.5)=1;
dat = filter(Num_LP, Den_LP, dat);   % Filter data between 0.5-10Hz
dat = filter(Num_HP, Den_HP, dat);

snleo = nlin_energy(dat, win*fs1);      % Smoothed Absolute Nonlinear Energy Operator
s1 = resample(snleo, fs2, fs1); s2 = s1;
s1(arx==1)=0;
epl = fs2*wc; 
for ii = epl+1:length(s1)                      % Baseline Correction
    r1 = ii-epl; r2 = ii-1; 
    s2(ii) = s1(ii)-min(s1(r1:r2));
   % s3(ii) = min(s1(r1:r2));
end
s2(s2<0)=0;
k = median(s2);
s2(s2>k) = k*(log(s2(s2>k))-(log(k)-1));

M = 200; N = length(s2);
th = quantile(s2, M);
erf1 = zeros(1, length(th)-1);
for ii = 1:length(th)-1
    dum = zeros(1, N);
    dum(s2>th(ii)) = 1;
    dum = process_ba_1ch(dum, fs2, 1, 2, 0);
    r1 = find(diff([0 dum 0]) == 1);
    erf1(ii) = length(r1);   % This is the cost function which threshold selection is based
end
q1 = find(erf1==max(erf1), 1, 'last');
erf1 = erf1(q1:end);
q2 = find([0 erf1]>[erf1 0]);
% indices where erf1 decreases
% A few checks and if not satisfied merely select the threshold as the
%median of the data
%The aim is to select the start of the longest period of constant (stable)
%ERF which is just burst number slected by a given threshold
if isempty(q2)==1
    th1 = th(M/2);
else
    th1 = th(q1+q2(find(diff(q2)==max(diff(q2)), 1, 'last'))-1);
    if isempty(th1)==1 
       th1 = th(M/2);
    end
end
dum = zeros(1, N);
dum(s2>th1) = 1;
bm = process_ba_1ch(dum, fs2, 1, 2, 0);

end

