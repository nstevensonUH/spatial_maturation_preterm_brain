function p = do_image(data, MyAtlas, flag)

braintype = 'real';
show_scalp = 0;       %Show scalp segment on headmodel
alphaval = 1;       %surface transparency

%==========================================================================  
% Load cx model
%==========================================================================  

% Load Atlas with areas and Circular xy (Nch 58)

% Load cx surface model       
load real_cx.mat;
load flat_cx.mat;
load scalp.mat;

Nvrt=length(real_cx.Vertices); % N of vertices
Np = length(MyAtlas.Parcels);  % N of parcels 

switch flag
    case 1 
        
        valr = 850/1000;
        set(gcf, 'Color', 'w');
        %data(isnan(data))=0;  %Replace possible nans with 0
        subplot('Position', [0.15 0.3 0.3*valr 0.3]); hold on;
        cdata = NaN*ones(Nvrt,1);
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

        %colormap('jet'); 
        ncolors = [32 32]
        CT0 = [1 0.765 0; 
    1 1 1];
 
% interpolate
x0 = linspace(0,1,size(CT0,1));
xf = linspace(0,1,ncolors(1));
CT1 = interp1(x0,CT0,xf,'linear');

CT0 = [1 1 1 ; 
    0.5 0 0.5];
% interpolate
x0 = linspace(0,1,size(CT0,1));
xf = linspace(0,1,ncolors(2));
CT2 = interp1(x0,CT0,xf,'linear');

CT = [repmat(CT1(1,:), 8,1) ; CT1(1:end-1,:) ; CT2];

% clamp values
CT = imclamp(CT); % your new color table

        % ncolors = [32 32];
        % CT0 = [1 0.765 0 ; 0.5 0 0.5];
        % %CT0 = 0.8*[51 255 255 ; 255 129 0]./255;
        % % interpolate
        % x0 = linspace(0,1,size(CT0,1));
        % xf = linspace(0,1,ncolors(1));
        % CT1 = interp1(x0,CT0,xf,'linear');
        % % CT0 = [1 1 1 ; 0.5 0 0.5];
        % % % interpolate
        % % x0 = linspace(0,1,size(CT0,1));
        % % xf = linspace(0,1,ncolors(2));
        % % CT2 = interp1(x0,CT0,xf,'linear');
        % CT = [repmat(CT1(1,:), 4,1) ; CT1 ; repmat(CT1(end,:), 4,1)];
        % % clamp values
        % CT = imclamp(CT); % your new color table
        colormap(CT)

        axis([-0.08 0.08 -0.065 0.055])
        axis off
        hold off    
        %Set top view
        view(90,0);
        zlim([-0.015 0.105])
        ylim([-0.06 0.06])
        %axis equal
        text(0, -0.06, 0.1, 'C', 'fontsize', 24, 'fontname', 'times')
        text(0.0, 0.04, 0, 'L', 'fontsize', 20, 'fontname', 'times')

        subplot('Position', [0.15 0 0.3*valr 0.3]); hold on;  
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

        subplot('Position', [0.45 0.05 0.5*valr 0.5]); hold on;
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
        
    case 2

        %==========================================================================  
        % Plot data on 3D brain surface
        %==========================================================================  
        valr = 500/1000;
        subplot('Position', [0.15 0.5 0.5*valr 0.5]); hold on;
        set(gcf, 'Color', 'w');

        %Set colours for the vertices according to the data of the parcels        
        cdata = NaN*ones(Nvrt,1);
        for n=1:length(data)    
            cdata(MyAtlas.Parcels{n})=repmat(data(n),length(MyAtlas.Parcels{n}),1);      
        end

        %Vertices without cdata:
        cdata(isnan(cdata)) = 0.25; %nanmin(data);    %Replace possible nans with min(data)
        %cdata(isnan(cdata)) = 0;              %Replace possible nans with zero

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
        colormap jet;    
        %axis([-0.08 0.08 -0.065 0.055])
        axis off
        hold off   
        %Set top view
        view(90,0); 
        zlim([-0.015 0.105])
        ylim([-0.06 0.06])
        %axis equal
        text(0.0, 0.04, 0, 'L', 'fontsize', 20, 'fontname', 'times')

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
        cscale = [0 0.5]; 
        caxis(cscale);
        %axis([-0.055 0.065 0 0.08])
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
        
    case 3
        
        cm = [0.9 0 0 ; 0 0.7 0];
        colormap(cm)
        valr = 850/1000;
        set(gcf, 'Color', 'w');
        
        subplot('Position', [0.125 0.05 0.425*valr 0.425]); hold on;
        cdata = NaN*ones(Nvrt,1);
        for n=1:length(data)    
            cdata(MyAtlas.Parcels{n})=repmat(data(n),length(MyAtlas.Parcels{n}),1);      
        end
        cdata(~isnan(cdata)) = 1; %nanmin(data);    %Replace possible nans with min(data)
        cdata(isnan(cdata)) = 0;              %Replace possible nans with zero

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

        
        axis off
        hold off    
        %Set top view
        view(-90,90);
        xlim([-0.06 0.06]);
        ylim([-0.06 0.06]);

%         zlim([-0.015 0.105])
%         ylim([-0.06 0.06])
%         %axis equal
        text(-0.06, -0.06, 0.1, 'B', 'fontsize', 24, 'fontname', 'times')
        text(0.055, 0.035, 0, 'FL', 'fontsize', 20, 'fontname', 'times')


        subplot('Position', [0.575 0.05 0.425*valr 0.425]); hold on;
        if strcmp(braintype,'real')>0
            % Show brain (real)
            p=patch('Vertices',real_cx.Vertices,'Faces',real_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none');        
        else
            % Show brain (flat)
            p=patch('Vertices',flat_cx.Vertices,'Faces',flat_cx.Faces,'CData',cdata,'CDataMapping','scaled','FaceColor','interp','FaceAlpha',alphaval,'EdgeColor','none'); 
        end

        
        p.SpecularColorReflectance = 0.9;
        p.SpecularExponent = 1;
        p.SpecularStrength = 0.2;
        p.DiffuseStrength = 0.5;
        p.AmbientStrength = 0.5;
        
        axis off
        hold off
        view(90,-90);
        xlim([-0.06 0.06]);
        ylim([-0.06 0.06]);

%         zlim([-0.015 0.105])
%         ylim([-0.06 0.06])
        
        light('Position',[ 0.0 -1.0  -0.04], 'Color', [0.4 0.4 0.4]);
        light('Position',[ 0.0  1.0  -0.04], 'Color', [0.4 0.4 0.4]);

        light('Position',[-1.0  0.0  -0.04], 'Color', [0.4 0.4 0.4]);
        light('Position',[ 1.0  0.0  -0.04], 'Color', [0.4 0.4 0.4]);

        light('Position',[ 0.0  0.0  0.0], 'Color', [0.4 0.4 0.4]);
        light('Position',[ 0.0  0.0  1.0], 'Color', [0.4 0.4 0.4]);

    case 4

        %load('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\burst_analysis\regional_analysis\src_n125\MyAtlas_n58.mat');
        %[~, t, ~] = xlsread('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JamesR\burst_analysis\regional_analysis\SA-axis vs EEG parcellation_V002.xlsx');
        %MyAtlas.SA = {t{7:64,3}}';
        
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
        
               
        [val, idx] = sort(dis, 'descend');
        nref = find(val<0,1);
        nval = [-(nref-1):1:-1 1:length(dis)-nref+1];
        nmap = [-12 -8   0 -11   0  -5  -9  -7 -13 -4  0 -10 -3  0 -6 -2 -1  0 3 6 7 5 0 4 2 0 1  8  9];
        nv = zeros(1,length(dum));
        nv(idx) = nmap; 
        nv(idx+29) = nmap; 
        
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
        
        % figure; set(gcf, 'Position', [800 200 1000 500])
        % set(gcf, 'Color', 'w');
        % %data(isnan(data))=0;  %Replace possible nans with 0
        %valr = 500/1000;
        subplot(2,3, 3); hold on 
        
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

        colormap(CT)
        
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

        h = colorbar;
        h.Label.String = 'SA (rank)';
        h.FontSize = 14;
        h.FontName = 'times';
        h.Ticks = [0 0.5*13/22 0.5];
        h.TickLabels = {'A', '0', 'S'};
        h.Label.FontSize = 16;
        text(0.055, 0.035, 0, 'FL', 'fontsize', 16, 'fontname', 'times')
        text(0.075, 0.035, 0, 'C', 'fontsize', 20, 'fontname', 'times')


        case 5 
        
        valr = 1;
        set(gcf, 'Color', 'w');
        data(data==0)=NaN;  %Replace possible nans with 0
        subplot('Position', [0.15 0.3 0.3*valr 0.3]); hold on;
        cdata = NaN*ones(Nvrt,1);
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

        %colormap('jet'); 
        ncolors = [32 32]
        CT0 = [1 0.765 0; 
                1 1 1];
 
        % interpolate
        x0 = linspace(0,1,size(CT0,1));
        xf = linspace(0,1,ncolors(1));
        CT1 = interp1(x0,CT0,xf,'linear');
        
        CT0 = [1 1 1 ; 
            0.5 0 0.5];
        % interpolate
        x0 = linspace(0,1,size(CT0,1));
        xf = linspace(0,1,ncolors(2));
        CT2 = interp1(x0,CT0,xf,'linear');
        
        CT = [repmat(CT1(1,:), 8,1) ; CT1(1:end-1,:) ; CT2];

        % clamp values
        CT = imclamp(CT); % your new color table

        % ncolors = [32 32];
        % CT0 = [1 0.765 0 ; 0.5 0 0.5];
        % %CT0 = 0.8*[51 255 255 ; 255 129 0]./255;
        % % interpolate
        % x0 = linspace(0,1,size(CT0,1));
        % xf = linspace(0,1,ncolors(1));
        % CT1 = interp1(x0,CT0,xf,'linear');
        % % CT0 = [1 1 1 ; 0.5 0 0.5];
        % % % interpolate
        % % x0 = linspace(0,1,size(CT0,1));
        % % xf = linspace(0,1,ncolors(2));
        % % CT2 = interp1(x0,CT0,xf,'linear');
        % CT = [repmat(CT1(1,:), 4,1) ; CT1 ; repmat(CT1(end,:), 4,1)];
        % % clamp values
        % CT = imclamp(CT); % your new color table
%        colormap(CT)

        axis([-0.08 0.08 -0.065 0.055])
        axis off
        hold off    
        %Set top view
        view(90,0);
        zlim([-0.015 0.105])
        ylim([-0.06 0.06])
        %axis equal
        text(0, -0.06, 0.1, 'C', 'fontsize', 24, 'fontname', 'times')
        text(0.0, 0.04, 0, 'L', 'fontsize', 20, 'fontname', 'times')

        subplot('Position', [0.15 0 0.3*valr 0.3]); hold on;  
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

        subplot('Position', [0.45 0.05 0.5*valr 0.5]); hold on;
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
        
    

        
end