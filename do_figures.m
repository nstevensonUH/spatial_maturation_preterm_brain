function do_figures()


load flist

% NEW ANALSYSIS SA AXIS ONLY, LATERALLY AVERAGED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load results_mse_v8_without
load MyAtlas_n58
[~, t, ~] = xlsread('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\burst_analysis\regional_analysis\SA-axis vs EEG parcellation_V002.xlsx');
MyAtlas.SA = {t{7:64,3}}';
idx = zeros(1,58);
idx([19 20 22 23 25 28 29 48 49 51 52 54 57 58]) = 2;  %sources at depth
val_ref = find(idx==0);

% NEW ANALSYSIS SA AXIS ONLY, FBA, 
dum = MyAtlas.Centroids;
dx = MyAtlas.SA;
for ii = 1:length(dx)
    if dx{ii} == 'S'
        b = contains(dx,'A');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = -mean(sqrt(sum((dum(rf,:) - pos).^2,2)));
        dval(ii) = 0;
    else
        b = contains(dx,'S');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = mean(sqrt(sum((dum(rf,:) - pos).^2,2)));
        dval(ii)=1;
    end
end
dis = dis(val_ref);
dval = dval(val_ref);
disr = MyAtlas.SArank(val_ref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% more or less correlated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r12 = corr(disr', MyAtlas.Centroids(val_ref,1));
out1 = cor_comp(r(2,1), r(2,4), r12, length(val_ref), 0.05); % signficantly more correlated
r12 = corr(abs(MyAtlas.Centroids(val_ref,2)), MyAtlas.Centroids(val_ref,1));
out2 = cor_comp(r(2,1), -r(2,2), r12, length(val_ref), 0.05); % not signficantly more correlated


pad = zeros(length(val_ref),1);
for ii = 1:length(val_ref)
    rf = 1:length(val_ref);
    rf(ii)=0; rf1 = find(rf>0);
    pad(ii) = median(err1(val_ref(ii),val_ref(rf1)));    
end

pav = zeros(length(val_ref),1);
for ii = 1:length(val_ref)
    rf = 1:length(val_ref);
    rf(ii)=0; rf1 = find(rf>0);
    pav(ii) = median(err2(val_ref(rf1),val_ref(ii)));    
end


cm(1,:) = [0.5 0 0.5]*0.9;
cm(2,:) = [1 0.765 0]*0.9;
figure; set(gcf, 'Position', [800 100 1000 850])
subplot(3,6,2:3); hold on;
x = dum(val_ref,1);
plot(100*x, pad, 'o', 'linewidth', 2, 'color', cm(1,:)); 
B = polyfit(x, pad', 1);
xx = linspace(min(x), max(x), 100);
plot(100*xx, polyval(B, xx), 'linewidth', 2, 'color', cm(2,:));
axis([-6 6 -2.1 1.5])
set(gca, 'fontsize', 16, 'position', [0.2 0.675 0.275 0.315], 'fontname', 'times', 'Xtick', [-5:2.5:5]); %, 'Xticklabel', [-0.05:0.025:0.05]*100)
ylabel('\Delta FBA (weeks)'); xlabel('frontal-occipital axis (cm)')
text(-5, 1.3, 'A', 'fontsize', 24, 'fontname', 'times'); grid off
[rd, pd] = corr(x, pad)
text(-1.25, -1.75, '{\it r} = 0.606; {\it p} < 0.001', 'fontsize', 14, 'fontname', 'times'); 


subplot(3,6,4:5); hold on;
x = abs(dum(val_ref,2));
plot(100*dum(val_ref,2), pad, 'o', 'linewidth', 2, 'color', cm(1,:)); 
B = polyfit(x, pad', 1);
xx1 = linspace(0, max(x), 100);
dmx = polyval(B, xx1);
plot(100*xx1, dmx, 'linewidth', 2, 'color', cm(2,:));
xx2 = linspace(-max(x), 0, 100);
plot(100*xx2, dmx(end:-1:1), 'linewidth', 2, 'color', cm(2,:));
axis([-5 5 -2.1 1.5])
set(gca, 'fontsize', 16, 'position', [0.6 0.675 0.275 0.315], 'fontname', 'times', 'Xtick', [-5:2.5:5]); %, 'Xticklabel', [-0.05:0.025:0.05]*100)
ylabel('\Delta FBA (weeks)'); xlabel('lateral axis (cm)')
text(-4, 1.3, 'B', 'fontsize', 24, 'fontname', 'times'); grid off
[rd, pd] = corr(x, pad)
text(-1.25, -1.5, '{\it r} = -0.487; {\it p} < 0.001', 'fontsize', 14, 'fontname', 'times'); 

age3D = zeros(1,length(dum));
age3D(val_ref)=pad;
rng_data = [min(age3D) max(age3D)];
age3D = (age3D-min(age3D));
age3D = age3D./max(age3D);
data = age3D/2; 

do_image(data, MyAtlas, 1);

cscale = [0 0.5]; 
caxis(cscale);
axis off
hold off   
view(270,90); 
xlim([-0.06 0.06]);
ylim([-0.07 0.05]);
h = colorbar;
h.Label.String = '\Delta FBA (weeks)';
h.FontSize = 14;
h.FontName = 'times';
h.Ticks = 0:0.1:0.5;
htl = round(linspace(-2, 1, length(h.Ticks))*10)/10;
h.TickLabels = num2str(htl');
h.Label.FontSize = 18;
text(0.055, 0.035, 0, 'FL', 'fontsize', 20, 'fontname', 'times')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% early vs late gradients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = zeros(1,58);
idx([19 20 22 23 25 28 29 48 49 51 52 54 57 58]) = 2;  %sources at depth
val_ref = find(idx==0);
dum = MyAtlas.Centroids;
dum = MyAtlas.Centroids;
dx = MyAtlas.SA;
for ii = 1:length(dx)
    if dx{ii} == 'S'
        b = contains(dx,'A');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = -mean(sqrt(sum((dum(rf,:) - pos).^2,2)));
        dval(ii) = 0;
    else
        b = contains(dx,'S');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = mean(sqrt(sum((dum(rf,:) - pos).^2,2)));
        dval(ii)=1;
    end
end
dis = dis(val_ref);
dval = dval(val_ref);
disr = MyAtlas.SArank(val_ref);

load results_mse_v8_without_vx.mat

errx1 = zeros(length(idx)); errx2 = errx1;
for z1 = 1:length(idx)
    for z2 = 1:length(idx)
        tr = out1{z1,z2};
        ts = out4{z1,z2};
        pm = out2{z1,z2};
        errx1(z1,z2) = mean(tr(pm<36)-ts(pm<36));
        errx2(z1,z2) = mean(tr(pm>=36)-ts(pm>=36));
    end
end


pad = zeros(length(val_ref),1); pad1 = pad; pad2 = pad;
for ii = 1:length(val_ref)
    rf = 1:length(val_ref);
    rf(ii)=0; rf1 = find(rf>0);
    pad(ii) = median(err1(val_ref(ii),val_ref(rf1)));    
    pad1(ii) = median(errx1(val_ref(ii),val_ref(rf1)));    
    pad2(ii) = median(errx2(val_ref(ii),val_ref(rf1)));    
end

[r, ~, c1, c2] = corrcoef(pad, dum(val_ref,1));
rs = [r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad, abs(dum(val_ref,2)));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad, dum(val_ref,3));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad, disr);
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];


[r, ~, c1, c2] = corrcoef(pad1, dum(val_ref,1));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad1, abs(dum(val_ref,2)));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad1, dum(val_ref,3));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad1, disr);
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];


[r, ~, c1, c2] = corrcoef(pad2, dum(val_ref,1));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad2, abs(dum(val_ref,2)));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad2, dum(val_ref,3));
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];
[r, ~, c1, c2] = corrcoef(pad2, disr);
rs = [rs ; r(1,2) c1(1,2) c2(1,2)];



figure; hold on;
x = dum(val_ref,1);
h1 = plot(100*x, pad, 'o', 'linewidth', 2); 
h1.MarkerFaceColor = h1.Color;
h2 = plot(100*x, pad1, 's', 'linewidth', 2); 
h2.MarkerFaceColor = h2.Color;
h3 = plot(100*x, pad2, '^', 'linewidth', 2); 
h3.MarkerFaceColor = h3.Color;
B1 = polyfit(x, pad', 1);
B2 = polyfit(x, pad1', 1);
B3 = polyfit(x, pad2', 1);
xx = linspace(min(x), max(x), 100);
plot(100*xx, polyval(B1, xx), 'linewidth', 2, 'color', h1.Color);
plot(100*xx, polyval(B2, xx), 'linewidth', 2, 'color', h2.Color);
plot(100*xx, polyval(B3, xx), 'linewidth', 2, 'color', h3.Color);
axis([-6 6 -2.1 1.5])
set(gca, 'fontsize', 16, 'fontname', 'times', 'Xtick', -5:2.5:5); %, 'Xticklabel', [-0.05:0.025:0.05]*100)
ylabel('\Delta FBA (weeks)'); xlabel('FO distance (cm)')
text(-5, 1.3, 'A', 'fontsize', 24, 'fontname', 'times'); grid on; box on;
legend([h1,h2,h3],'all','< 36 w PMA','\geq 36 w PMA')
[rd, pd] = corr(x, pad)
%text(-1.25, -1.75, '{\it r} = 0.606; {\it p} < 0.001', 'fontsize', 14, 'fontname', 'times'); 

%
% FIGURE 3
%

cm(1,:) = [0.5 0 0.5]*0.9;
cm(2,:) = [1 0.765 0]*0.9;
figure; set(gcf, 'Position', [800 100 800 400])
subplot(2,3, [1 4]); hold on;
x = disr;
plot(x, pad, 'o', 'linewidth', 2, 'color', cm(1,:)); hold on;
B = polyfit(x, pad', 1);
xx = linspace(min(x), max(x), 100);
plot(xx, polyval(B, xx), 'linewidth', 2, 'color', cm(2,:));
axis([-10 15 -2.1 1.5])
set(gca, 'fontsize', 16, 'position', [0.10 0.155 0.3 0.8], 'fontname', 'times', 'Xtick', [-10 0 15], 'Xticklabel', {'S','0','A'})
ylabel('\Delta FBA (weeks)'); xlabel('SA rank')
text(-8, 1.3, 'A', 'fontsize', 24, 'fontname', 'times'); grid off

do_image([], MyAtlas, 4)
set(gca, 'position', [0.725 0.4 0.25 0.5])


subplot(2,3, [2 5]); hold on;
pv = ranksum(pad(dval==1), pad(dval==0));
plot_box_violin(0, pad(dval==0), cm(2,:), [1 1], -1); % Sensori-motor
plot_box_violin(1, pad(dval==1), cm(1,:), [1 1], 1); % Association
ylabel('\Delta FBA (weeks)'); 
set(gca, 'fontsize', 16, 'position', [0.5 0.155 0.2 0.8], 'fontname', 'times', 'Xtick', [0 1], 'Xticklabel', {'S', 'A'})
axis([-1 2 -2.1 1.5])
plot([0 1], [1.2 1.2], 'k', 'linewidth', 2)
text(0.45, 1.3, '*', 'fontsize', 24)
text(-1, 1.3, 'B', 'fontsize', 24, 'fontname', 'times'); grid off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIGURES FOR SUPPLEMENTAL MATERIAL
%   S1: regions ignored for analysis    
%   S2: ranking of SA axis
%   S3: differences in gradient between FO and SA axis
%   S4: gradients for FO axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% S1: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load MyAtlas_n58
load real_cx.mat;
load flat_cx.mat;
load scalp.mat;
braintype = 'real';     
Nvrt=length(real_cx.Vertices); % N of vertices
Np = length(MyAtlas.Parcels);  % N of parcels 
show_scalp = 0;       %Show scalp segment on headmodel
alphaval = 1;       %surface transparency

idx = zeros(1,58);
idx([19 20 22 23 25 28 29 48 49 51 52 54 57 58]) = 2;  %sources at 
dref = idx; dref = find(idx==2);
val_ref = find(idx==0);

pas = mean(err4);

figure; set(gcf, 'Position', [800 100 1000 850])
subplot(3,6,3:4); hold on;
plot(pas, 'color', [0 0 0.5]); 
p1 = plot(dref, pas(dref), 'o')
axis([1 58 1.2 2.6])
set(gca, 'fontsize', 16, 'position', [0.375 0.675 0.275 0.315], 'fontname', 'times')
ylabel('Systematic Error'); xlabel('Parcel')
text(5, 2.2, 'A', 'fontsize', 24, 'fontname', 'times'); grid off
legend(p1, 'Excluded')

data = zeros(Np,1);
data(idx==2) = NaN;

do_image(data, MyAtlas, 3);

light('Position',[ 0.0 -1.0  -0.04], 'Color', [0.4 0.4 0.4]);
light('Position',[ 0.0  1.0  -0.04], 'Color', [0.4 0.4 0.4]);

h = colorbar;
h.FontSize = 14;
h.Ticks = ([0 0.5]);
h.TickLabels = {'invalid','valid'};
h.FontSize = 18;
h.FontName = 'times';
text(0.055, 0.035, 0, 'FL', 'fontsize', 20, 'fontname', 'times')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  S2: ranking of SA axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\burst_analysis\regional_analysis\src_n125\MyAtlas_n58.mat');
% [~, t, ~] = xlsread('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\burst_analysis\regional_analysis\SA-axis vs EEG parcellation_V002.xlsx');
% MyAtlas.SA = {t{7:64,3}}';

dum = MyAtlas.Centroids;
dx = MyAtlas.SA; sa = zeros(1,length(dx)/2); dis = zeros(1, length(dx)/2);
for ii = 1:length(dx)/2
    if dx{ii} == 'S'
        b = contains(dx,'A');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = -(mean(sqrt(sum((dum(rf,:) - pos).^2,2))));
        sa(ii) = 0;
    else
        b = contains(dx,'S');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = mean(sqrt(sum((dum(rf,:) - pos).^2,2)));
        sa(ii) = 1;
    end
end

% initial faffing line line up with Syndor et al. style parcellation
% figure;
% idx = zeros(1,58);
% idx([19 20 22 23 25 28 29 48 49 51 52 54 57 58]) = 2;  %sources at 
% val_ref = find(idx==0); sax = sa(val_ref(val_ref<29)); vr = val_ref(val_ref<29);
% dm = dum(val_ref,:);
% figure; plot(-dm(sax==0,1), dm(sax==0,3), 'o', 'linewidth', 2)
% hold on; plot(-dm(sax==1,1), dm(sax==1,3), '+', 'linewidth', 2)
% plot(-dum(idx==2,1), dum(idx==2,3), 'k^', 'linewidth', 2)
% set(gca, 'fontsize', 14)
% legend('S','A')
% text(-0.05, 0.005, 'front', 'fontsize', 14)
% axis([-0.06 0.06 0 0.08])
% grid on
% xlabel('frontal-occpital axis')
% ylabel('inferior-superior axis')


[val, idx] = sort(dis, 'descend');
nref = find(val<0,1);
nval = [-(nref-1):1:-1 1:length(dis)-nref+1];
%figure; subplot(1,2,1) 
% for ii = 1:length(dum)/2
%   %  if ~isempty(find(vr==idx(ii)))
%     text(-dum(idx(ii),1), dum(idx(ii),3), num2str(nval(ii)), 'fontsize', 14)
%    % end
% end
% axis([-0.06 0.06 0 0.08])
% set(gca, 'fontsize', 14)
% text(-0.05, 0.005, 'front', 'fontsize', 14)
% axis([-0.06 0.06 0 0.08])
% grid on
% xlabel('frontal-occpital axis')
% ylabel('inferior-superior axis')
% title('old ranking')

%nmap = [-18 -14 -15 -17 -11 -4 -13 -8 -16 -3 -9 -12 -2 -10 -5 -7 -6 -1 5 9 8 7 1 6 4 3 2 10 11];
%       [-18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 1 2 3 4 5 6 7 8 9 10 11] 
nmap = [-12 -8   0 -11   0  -5  -9  -7 -13 -4  0 -10 -3  0 -6 -2 -1  0 3 6 7 5 0 4 2 0 1  8  9];
% subplot(1,2,2)
% for ii = 1:length(dum)/2
%     if ~isempty(find(vr==idx(ii)))
%     text(-dum(idx(ii),1), dum(idx(ii),3), num2str(nmap(ii)))
%     end
% end
% axis([-0.06 0.06 0 0.08])
% set(gca, 'fontsize', 14)
% text(-0.05, 0.005, 'front', 'fontsize', 14)
% axis([-0.06 0.06 0 0.08])
% grid on
% xlabel('frontal-occpital axis')
% ylabel('inferior-superior axis')
% title('new ranking')

nv = zeros(1,length(dum));
nv(idx) = nmap; 
nv(idx+29) = nmap; 

% figure; 
% for ii = 1:length(dum)/2
%     text(-dum(ii,1), dum(ii,3), num2str(nv(ii)))
% end
% axis([-0.06 0.06 0 0.08])

dis = nv;
rng_data = [min(dis) max(dis)]
dis = (dis-min(dis));
dis = dis./max(dis);
data = dis/2;

braintype = 'real'
alphaval = 1;       %surface transparency
% Load cx surface model       
load real_cx.mat;
load flat_cx.mat;
load scalp.mat;   
Nvrt=length(real_cx.Vertices); % N of vertices
Np = length(MyAtlas.Parcels);  % N of parcels 


ncolors = [14 10];
% some approximate breakpoint colors
% purple-white-yellow
CT0 = [0.5 0 0.5; 
    1 1 1];
% interpolate
x0 = linspace(0,1,size(CT0,1));
xf = linspace(0,1,ncolors(1));
CT1 = interp1(x0,CT0,xf,'linear');

CT0 = [1 1 1 ; 
    1 0.765 0];
% interpolate
x0 = linspace(0,1,size(CT0,1));
xf = linspace(0,1,ncolors(2));
CT2 = interp1(x0,CT0,xf,'linear');

CT = [CT1(1:end-1,:) ; CT2];

% clamp values
CT = imclamp(CT); % your new color table
% display the CT as an image for sake of visualization
%image(permute(flipud(CT),[1 3 2]))

figure; set(gcf, 'Position', [800 200 1000 500])
set(gcf, 'Color', 'w');
%data(isnan(data))=0;  %Replace possible nans with 0
valr = 500/1000;
subplot('Position', [0.15 0.5 0.5*valr 0.5]); hold on;
cdata = zeros(Nvrt,1);
for n=1:length(data)    
    cdata(MyAtlas.Parcels{n})=repmat(data(n),length(MyAtlas.Parcels{n}),1);      
end
if strcmp(braintype,'real')>0
    % Show brain (real)
    p=patch('Vertices',real_cx.Vertices,'Faces',real_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none');        
else
    % Show brain (flat)
    p=patch('Vertices',flat_cx.Vertices,'Faces',flat_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none'); 
end

light('Position',[ 0.0 -1.0  0.04], 'Color', [0.4 0.4 0.4]);
light('Position',[ 0.0  1.0  0.04], 'Color', [0.4 0.4 0.4]);

light('Position',[-1.0  0.0  0.04], 'Color', [0.4 0.4 0.4]);
light('Position',[ 1.0  0.0  0.04], 'Color', [0.4 0.4 0.4]);

light('Position',[ 0.0  0.0  0.0], 'Color', [0.4 0.4 0.4]);
light('Position',[ 0.0  0.0  1.0], 'Color', [0.4 0.4 0.4]);

p.SpecularColorReflectance = 0.9;
p.SpecularExponent = 1;
p.SpecularStrength = 0.2;
p.DiffuseStrength = 0.5;
p.AmbientStrength = 0.5;

%colormap('winter');    
%colormap('spring');    
%%%%%%

axis([-0.08 0.08 -0.065 0.055])
axis off
hold off    
%Set top view
view(90,0);
zlim([-0.015 0.105])
ylim([-0.06 0.06])
%axis equal
%text(0, -0.06, 0.1, 'B', 'fontsize', 24, 'fontname', 'times')
text(0.0, 0.04, 0, 'L', 'fontsize', 20, 'fontname', 'times')

colormap(CT)

subplot('Position', [0.15 0 0.5*valr 0.5])  
if strcmp(braintype,'real')>0
    % Show brain (real)
    p=patch('Vertices',real_cx.Vertices,'Faces',real_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none');        
else
    % Show brain (flat)
    p=patch('Vertices',flat_cx.Vertices,'Faces',flat_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none'); 
end
light('Position',[ 0.0 -1.0  0.04], 'Color', [0.4 0.4 0.4]);
light('Position',[ 0.0  1.0  0.04], 'Color', [0.4 0.4 0.4]);

light('Position',[-1.0  0.0  0.04], 'Color', [0.4 0.4 0.4]);
light('Position',[ 1.0  0.0  0.04], 'Color', [0.4 0.4 0.4]);

light('Position',[ 0.0  0.0  0.0], 'Color', [0.4 0.4 0.4]);
light('Position',[ 0.0  0.0  1.0], 'Color', [0.4 0.4 0.4]);

p.SpecularColorReflectance = 0.9;
p.SpecularExponent = 1;
p.SpecularStrength = 0.2;
p.DiffuseStrength = 0.5;
p.AmbientStrength = 0.5;

%==========================================================================  
% Set color map and scale (min max)
%==========================================================================  
cscale = [0 0.5]; 
caxis(cscale);
axis([-0.055 0.065 0 0.08])
axis off
hold off   
%Set top view
view(180,0); 
zlim([-0.025 0.095])
xlim([-0.06 0.06])
%axis equal
text(0.06, 0.00, 0.01, 'F', 'fontsize', 20, 'fontname', 'times')
subplot('Position', [0.5 0.1625 0.625*valr 0.625]); hold on;
if strcmp(braintype,'real')>0
    % Show brain (real)
    p=patch('Vertices',real_cx.Vertices,'Faces',real_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none');        
else
    % Show brain (flat)
    p=patch('Vertices',flat_cx.Vertices,'Faces',flat_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none'); 
end

light('Position',[ 0.0 -1.0  0.04], 'Color', [0.4 0.4 0.4]);
light('Position',[ 0.0  1.0  0.04], 'Color', [0.4 0.4 0.4]);

light('Position',[-1.0  0.0  0.04], 'Color', [0.4 0.4 0.4]);
light('Position',[ 1.0  0.0  0.04], 'Color', [0.4 0.4 0.4]);

light('Position',[ 0.0  0.0  0.0], 'Color', [0.4 0.4 0.4]);
light('Position',[ 0.0  0.0  1.0], 'Color', [0.4 0.4 0.4]);

p.SpecularColorReflectance = 0.9;
p.SpecularExponent = 1;
p.SpecularStrength = 0.2;
p.DiffuseStrength = 0.5;
p.AmbientStrength = 0.5;

%==========================================================================  
% Set color map and scale (min max)
%==========================================================================  
cscale = [0 0.5]; 
caxis(cscale);
axis([ -0.08 0.08 -0.07 0.05])
axis off
hold off   
%set(gca, 'Position', [0.5 0 0.5 1])  
%Set top view
view(270,90); 
xlim([-0.06 0.06]);
ylim([-0.07 0.05]);
%axis equal 
h = colorbar;
h.Label.String = 'SA (rank)';
h.FontSize = 14;
h.FontName = 'times';
h.Ticks = [0 0.5*13/22 0.5];
h.TickLabels = {'A', '0', 'S'};
h.Label.FontSize = 18;
text(0.055, 0.035, 0, 'FL', 'fontsize', 20, 'fontname', 'times')

%
%
% S4: all axis
%
%

% idx = zeros(1,58);
% idx([19 20 22 23 25 28 29 48 49 51 52 54 57 58]) = 2;  %sources at depth
% val_ref = find(idx==0);
% 
% % NEW ANALSYSIS SA AXIS ONLY, FBA, 
% dum = MyAtlas.Centroids;
% pad = zeros(length(val_ref),1);
% for ii = 1:length(val_ref)
%     rf = 1:length(val_ref);
%     rf(ii)=0; rf1 = find(rf>0);
%     pad(ii) = median(err1(val_ref(ii),val_ref(rf1)));    
% end
% figure; set(gcf, 'position', [200 400 1380 460])
% subplot(1,3,1); hold on;
% plot(dum(val_ref,1)*100, pad, 'o', 'linewidth', 2); hold on;
% B = polyfit(dum(val_ref,1), pad', 1);
% x = linspace(-6, 6, 100)/100;
% plot(x*100, polyval(B, x), 'linewidth', 2);
% axis([-6 6 -2 1.5])
% set(gca, 'fontsize', 16, 'fontname', 'times', 'Xtick', [-6:2:6], 'position', [0.06 0.135 0.275 0.8])
% ylabel('\Delta FBA (weeks)'); xlabel('Frontal-Occipital Distance (cm)'); grid on
% text(-5.5, 1.25, 'A', 'fontname', 'times', 'fontsize', 24)
% 
% subplot(1,3,2); hold on;
% plot(dum(val_ref,2)*100, pad, 'o', 'linewidth', 2); hold on;
% B = polyfit(dum(val_ref,2), pad', 1);
% x = linspace(-4, 4, 100)/100;
% plot(x*100, polyval(B, x), 'k', 'linewidth', 2);
% axis([-4.5 4.5 -2 1.5])
% set(gca, 'fontsize', 16, 'fontname', 'times', 'Xtick', [-4:2:4], 'position', [0.4 0.135 0.275 0.8])
% ylabel('\Delta FBA (weeks)'); xlabel('Lateral Distance (cm)'); grid on;
% text(-4, 1.25, 'B', 'fontname', 'times', 'fontsize', 24)
% 
% subplot(1,3,3); hold on;
% plot(dum(:,3)*100, pad, 'o', 'linewidth', 2); hold on;
% B = polyfit(dum(:,3), pad', 3);
% x = linspace(0.5, 8, 100)/100;
% plot(x*100, polyval(B, x), 'linewidth', 2);
% axis([0.5 8 -0.9 1.6])
% set(gca, 'fontsize', 16,  'fontname', 'times')
% xlabel('Inferior-Superior Distance (cm)'); grid on;
% 
% %
% %
% % S4: 
% %
% %
% 


load results_mse_v8_without
load all_data

idx = zeros(1,58);
idx([19 20 22 23 25 28 29 48 49 51 52 54 57 58]) = 2;  %sources at depthval_ref = find(idx==0);
val_ref = find(idx==0);
dum = MyAtlas.Centroids;
dx = MyAtlas.SA;
for ii = 1:length(dx)
    if dx{ii} == 'S'
        b = contains(dx,'A');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = -mean(sqrt(sum((dum(rf,:) - pos).^2,2)));
    else
        b = contains(dx,'S');
        rf = find(b==1);
        pos = dum(ii,:);
        dis(ii) = mean(sqrt(sum((dum(rf,:) - pos).^2,2)));
    end
end
dum  = dum(val_ref,:);
dis = dis(val_ref);
dis = round(dis*10000)/10000;
[~, idx] = sort(dis);
ud = unique(dis);
disr = zeros(size(dis));
for ii = 1:length(ud)
    aa = find(dis==ud(ii));
    disr(aa) = find(idx==aa(1));
end

% NTOE TO SELF USE REFERENCES FROM OUT4 for each file rather than an
% rmoutliers for fairer comparison
% Gradient in Burst Sharpness
pma = pma2';
val = zeros(20, 5);
for tt = 1:length(fref)
    ii = fref(tt);
    Z1 = size(fv2); Z2 = size(fv2{1});
    fts = zeros(Z1(2), Z2(1));
    for z1 = 1:Z1(2)
        fts(z1,:) = fv2{z1}(:,ii)';
    end
    % INITIAL FILTER STAGE
    erx = zeros(length(val_ref)); cx = zeros(1, length(val_ref));
    for ii = 1:length(val_ref)
        for jj = 1:length(val_ref)
            erx(ii,jj) = mean(fts(out4{val_ref(ii), val_ref(jj)},val_ref(ii))-fts(out4{val_ref(ii), val_ref(jj)},val_ref(jj)));
        end
        cx(ii) = corr(pma(out4{val_ref(ii), val_ref(ii)}), fts(out4{val_ref(ii), val_ref(ii)}, val_ref(ii)))';
    end

    %[r1x, ~, c1x, c2x] = corrcoef(dum(val_ref,1), mean(erx));
    val(tt,:) = [tt median(cx) corr(dum(:,1), mean(erx)') corr(abs(dum(:,2)), mean(erx)') corr(disr', mean(erx)')] ;
end

M1 = median(median(cellfun(@length, out4)));
M2 = 44;
x1 = linspace(0,1,M1);
x2 = linspace(0,1,M2);
r1x = zeros(1,200); r2x = r1x; p1x = r1x; p2x = r1x; c1 = 1;
for ii = 0.5:0.25:50
    r1 = zeros(1,4000); r2 = r1; p1 = r1; p2 = r1;
    for jj = 1:4000
        [r1(jj), p1(jj)] = corr(x1', (x1+(randn(1,M1))*sqrt(ii/10))') ;
        [r2(jj), p2(jj)] = corr(x2', (x2+(randn(1,M2))*sqrt(ii/10))') ;
    end
    r1x(c1) = mean(r1); p1x(c1) = mean(p1);
    r2x(c1) = mean(r2); p2x(c1) = mean(p2);
    c1 = c1+1;
end

rx = [1:100 101:4:200];
beta1 = nlinfit(r1x(rx), p1x(rx), @gfit, 0.01);
rx = [1:50 51:8:200];
beta2 = nlinfit(r2x(rx), p2x(rx), @gfit, 0.07);
nr = linspace(0,1,10000);
th1 = nr(find(gfit(beta1, nr)<=0.05/20, 1));
th2 = nr(find(gfit(beta2, nr)<=0.05/20, 1));

th1x = nr(find(gfit(beta1, nr)<=0.05, 1));
th2x = nr(find(gfit(beta2, nr)<=0.05, 1));

[d1, d2] = sort(abs(val(:,3)), 'descend');
rv1 = r(2,1); 
rv2 = median(diag(err3(val_ref,val_ref)));

cm(1,:) = [100, 149, 237]/255;
cm(2,:) = [128, 0, 0]/255;

cm1(1,:) = [100, 149, 237]/255;
cm1(2,:) = [129, 0, 129]/255;
cm1(3,:) = [0, 196, 0]/255;
cm1(4,:) = [255  196 0]/255;

cdx = [2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 4 3 3];

figure; set(gcf, 'Position', [700 550 900 420]); subplot(2,1,1); hold on;
for ii = 1:20
    plot(ii, d1(ii), 'o', 'markersize', 14, 'MarkerEdgeColor', cm1(cdx(d2(ii)),:), 'MarkerFaceColor', cm1(cdx(d2(ii)),:))
end
plot(-1, rv1, 's', 'markersize', 14, 'MarkerEdgeColor', cm(2,:), 'MarkerFaceColor', cm(2,:)); 
plot([0 20.5], [th2 th2], 'k--')
h1 = plot([-2 0], [th1 th1], 'k--')
axis([-2 21 0 1])
set(gca, 'Fontsize', 16, 'Fontname', 'times', 'Ytick', [0:0.2:1], 'Xtick', [-1 1:20], 'Xticklabel', ['FBA' ; num2cell(d2)])
set(gca, 'Position', [0.1 0.595 0.85 0.37])
ylabel('|{\it r} | (distance)')
text(0, 0.8, 'A', 'Fontsize', 24, 'Fontname', 'times')
legend(h1, 'p=0.05')
subplot(2,1,2); hold on;
for ii = 1:20
    plot(ii, abs(val(d2(ii),2)), 'o',  'markersize', 14, 'MarkerEdgeColor', cm1(cdx(d2(ii)),:), 'MarkerFaceColor', cm1(cdx(d2(ii)),:))
end
%plot(abs(val(d2,2)), 'o', 'MarkerEdgeColor', cm(1,:), 'MarkerFaceColor', cm(1,:));
plot(-1, rv2, 's', 'markersize', 14, 'MarkerEdgeColor', cm(2,:), 'MarkerFaceColor', cm(2,:));
h1 = plot([0 20.5], [th2x th2x], 'k--')
plot([-2 0], [th1x th1x], 'k--')
set(gca, 'Fontsize', 16, 'Fontname', 'times', 'Ytick', [0:0.2:1], 'Xtick', [-1 1:20], 'Xticklabel', ['FBA' ; num2cell(d2)])
ylabel('|{\it r} | (age)')
axis([-2 21 0 1])
set(gca, 'Position', [0.1 0.15 0.85 0.37])
xlabel('Feature of Cortical Activity ID')
text(0, 0.8, 'B', 'Fontsize', 24, 'Fontname', 'times')


h2 = plot(-1,-1, 'o', 'markersize', 14, 'MarkerEdgeColor', cm1(1,:), 'MarkerFaceColor', cm1(1,:));
h3 = plot(-1,-1, 'o', 'markersize', 14, 'MarkerEdgeColor', cm1(2,:), 'MarkerFaceColor', cm1(2,:));
h4 = plot(-1,-1, 'o', 'markersize', 14, 'MarkerEdgeColor', cm1(3,:), 'MarkerFaceColor', cm1(3,:));
h5 = plot(-1,-1, 'o', 'markersize', 14, 'MarkerEdgeColor', cm1(4,:), 'MarkerFaceColor', cm1(4,:));
h6 =  plot(-1,-1, 's', 'markersize', 14, 'MarkerEdgeColor', cm(2,:), 'MarkerFaceColor', cm(2,:));
legend([h1, h6, h2, h3, h4, h5], 'p=0.05', 'FBA', 'Burst', 'Amplitude', 'Spectral', 'Information')


flist{1} = 'Median Amplitude Envelope'; %Amplitude Envelope is estimate using the analytic associate of a signal (via Hilbert transform)
flist{2} = '5th Percentile Amplitude Envelope';
flist{3} = '25th Percentile Amplitude Envelope';
flist{4} = '75th Percentile Amplitude Envelope';
flist{5} = '95th Percentile Amplitude Envelope';
flist{6} = '5th Percentile Inter-burst Interval';
flist{7} = '50th Percentile Inter-burst Interval';
flist{8} = '95th Percentile Inter-burst Interval';
flist{9} = '5th Percentile Burst Duration';
flist{10} = '50th Percentile Burst Duration';
flist{11} = '95th Percentile Burst Duration';
flist{12} = 'Burst Skewness/Symmetry (15.625-125ms)';
flist{13} = 'Burst Skewness/Symmetry (125-250ms)';
flist{14} = 'Burst Skewness/Symmetry (250-500ms)';
flist{15} = 'Burst Skewness/Symmetry (0.5-1s)';
flist{16} = 'Burst Skewness/Symmetry (1-2s)';
flist{17} = 'Burst Skewness/Symmetry (2-4s)';
flist{18} = 'Burst Skewness/Symmetry (15.625ms-128s)';
flist{19} = 'Burst Kurtosis/Sharpness (15.625-125ms)';
flist{20} = 'Burst Kurtosis/Sharpness (125-250ms)';
flist{21} = 'Burst Kurtosis/Sharpness (250-500ms)';
flist{22} = 'Burst Kurtosis/Sharpness (0.5-1s)';
flist{23} = 'Burst Kurtosis/Sharpness (1-2s)';
flist{24} = 'Burst Kurtosis/Sharpness (2-4s)';
flist{25} = 'Burst Kurtosis/Sharpness (15.625ms-128s)';
flist{26} = 'Slope Burst Duration vs Burst Area';
flist{27} = 'Intercept Burst Duration vs Burst Area';
flist{28} = 'Slope Burst Duration vs Skewness';
flist{29} = 'Intercept Burst Duration vs Skewness';
flist{30} = 'Slope Burst Duration vs Kurtosis';
flist{31} = 'Intercept Burst Duration vs Kurtosis';
flist{32} = 'Suppression Curve';
flist{33} = 'Mean Burst Duration';
flist{34} = 'Standard Deviation Burst Durations';
flist{35} = 'Power Law of Burst Duration Distribution alpha 1';
flist{36} = 'Power Law of Burst Duration Distribution alpha 2';
flist{37} = 'Power Law of Burst Duration Distribution alpha 3';
flist{38} = 'Power Law of Burst Area Distribution alpha 1';
flist{39} = 'Power Law of Burst Area Distribution alpha 2';
flist{40} = 'Power Law of Burst Area Distribution alpha 3';
flist{41} = 'Sample Entropy';
flist{42} = 'Log PSD Energy (wideband)';
flist{43} = 'Relative Delta 1 power';
flist{44} = 'Relative Delta 2 power';
flist{45} = 'Relative Theta power';
flist{46} = 'Relative Alpha power';
flist{47} = 'Relative Beta power';
flist{48} = 'Spectral Slope';
flist{49} = 'Spectral Edge Frequency (90th Percentile)';
