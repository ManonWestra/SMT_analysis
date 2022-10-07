%% extract_confined
% confinement analysis of trajectories longer than minSteps
% circles with radius of confinement zone are drawn in black (immobile) or
% red (mobile)
% INPUT:
%   trajectories.mat
% OUTPUT:
%   confinement.mat
%       confinement_data:     confinement data for all tracks            
%       confined_data:        confinement data for confined tracks
%       confinement_zones:    data of all confinement zones
%       fraction_confined:    fraction of tracks with confinement
%   updated trajectory struct with L, Dr, R, Tdwell, Dmax, Lc per track
%   confinement heat maps
%
% Manon Westra, 2022

clear;clc;
load trajectories.mat;

extract_confined_param  % load parameters

% ANALYSIS
index_tracklength = find([trajectory.length] >= minimum);

% which trajectories to analyze?
if strcmp(selected_tracks,'mobile')
    a = mobile;
    b = immobile;
elseif strcmp(selected_tracks,'immobile')
    a = immobile;
    b = mobile; 
elseif strcmp(selected_tracks,'all')
    a = index_tracklength;
else
    error('Error: Define which tracks to analyze.')
end
  
% confinement analysis of trajectories with length minSteps or longer
CI = [];                        % confinement index
MeanSqD = [];                   % msd first [minimum] time lags in um2

confinement_data = zeros(length(a),14);
L_all = [];
cz_index = [];
cz_x = [];
cz_y = [];
cz_rad = [];
cz_dwell = [];
cz_trans = logical([]);
Di = [];

for tj = 1:length(a)
    clear Do mDo
    xy = trajectory(a(tj)).xy;              % in um
    fr = trajectory(a(tj)).frame;           % framenumbers in track
    Dins = trajectory(a(tj)).Dinst;    % Dinst um2/s
    A = [xy fr];
    [L, Dr, Dmax] = ConfinementIndex(A, Sm, dt, length(xy),npointsMSDset, Smin, Dins); % calculate running L (confinement index) and D in um2/s
    [Res, Tdwell, T_Lc] = RfromL(A, L, alfa, Lcm, Tc, Lmax, Lcmax); % Res: (1)radius (2) x (3) y (4) dwell time (5) start (6) end
        
    cregions = size(Res);
    Di_av = [];
    
    confinement_data(tj, 1) = a(tj);                % trajectory number
    confinement_data(tj, 2) = length(xy);           % trajectory length  
    confinement_data(tj, 3) = T_Lc(1);              % total time
    
    if cregions(1) >0  
        vector_confined = zeros(length(Dr),1);
        Do = NaN;
        nDo = 0;
        for c = 1:cregions(1)
            if c == 1
                if Res(c,5) >= np  % if part before confinement is equal or bigger than np points (to calculate Dinst)
                    nDo = nDo + 1;
                    msd = MSD_new(A(1:Res(c,5),:),dt);
                    checkD = Dinst(msd(:,4),msd(:,1),np);
                    if isempty(checkD)
                        Do(nDo) = NaN;
                    else
                        Do(nDo) = checkD;
                    end
                end
            else
                if Res(c,5)-Res(c-1,6) >= np
                    nDo = nDo + 1;
                    msd = MSD_new(A(Res(c-1,6):Res(c,5),:),dt); % if part between confinement is equal or bigger than np points
                    checkD = Dinst(msd(:,4),msd(:,1),np);
                    if isempty(checkD)
                        Do(nDo) = NaN;
                    else
                        Do(nDo) = checkD;
                    end
                end
            end
            vector_confined(Res(c,5):Res(c,6))=1;
            msd = MSD_new(A(Res(c,5):Res(c,6),:),dt);
            Di = [Di; Dinst(msd(:,4),msd(:,1),np)];  % Dinst inside confinement regions for all confinement zones
            Di_av = [Di_av; Dinst(msd(:,4),msd(:,1),np)];  % Dinst inside confinement regions per track
        end
        if size(L,1)-Res(c,6) >= np   % if part from last confinement to end track is equal or bigger than np points
            nDo = nDo +1;
            msd = MSD_new(A(Res(c,6):size(L,1),:),dt);
            checkD = Dinst(msd(:,4),msd(:,1),np);
            if isempty(checkD)
                Do(nDo) = NaN;
            else
                Do(nDo) = checkD;
            end
        end
        if ~isempty(Do)
            mDo = median(Do,'omitnan');
            confinement_data(tj, 9) = mDo;          % median Dinst outside confinement regions um2/s
            confinement_data(tj, 12) = log10(mDo);  % log10 um2/s
        end        
        vector_confined = logical(vector_confined);
        Din = mean(Dr(vector_confined));
        Dout = mean(Dr(~vector_confined));
        
        confinement_data(tj, 4)  = Tdwell(1);               % total dwell time confined
        confinement_data(tj, 6)  = cregions(1);             % number of confinement regions
        confinement_data(tj, 7)  = mean(Res(:, 1)).*1000;   % average radius confinement regions (nm)
        confinement_data(tj, 8)  = median(Di_av);           % median Dinst inside confinement regions um2/s
        confinement_data(tj, 10) = Tdwell(2);               % percentage confined of total time
        confinement_data(tj, 11) = log10(median(Di_av));    % log10 um2/s
        confinement_data(tj, 13) = Dmax;                    % max D used for confinement index
        confinement_data(tj, 14) = T_Lc(2);                 % threshold L for this track
        Res(:,1:3)=Res(:,1:3)*1000;                         % radius, x and y from um to nm
        for c = 1:cregions(1)
            cz_index = [cz_index; a(tj)];   % make indexlist with confined tracks
            cz_x = [cz_x; Res(c,2)];        % list of x center confinement zones (nm)
            cz_y = [cz_y; Res(c,3)];        % list of y center confinement zones (nm)
            cz_rad = [cz_rad; Res(c,1)];    % list of confinement radii (nm)
            cz_dwell = [cz_dwell; Res(c,4)]; % list of dwell time per confinement zone (s)
            if Res(c,5) > 1 && Res(c,6) < size(L,1)
                cz_transtemp = true;
            else
                cz_transtemp = false;
            end
            cz_trans = [cz_trans; cz_transtemp]; % list with 1 if a transient cz (not begin or end)
        end
    end
    confinement_data(tj, 5) = confinement_data(tj, 3)-confinement_data(tj, 4);  % total time free

    trajectory(a(tj)).L = L(:,2);                                   % confinement index over time
    trajectory(a(tj)).Dr = Dr;                                      % diffusion coefficient um2/s over time
    trajectory(a(tj)).R = Res;                                      % results confinement in nm
    trajectory(a(tj)).Tdwell = Tdwell(1);                           % total dwell time confined
    trajectory(a(tj)).Dmax = Dmax;                                  % max D used for confinement calculation
    trajectory(a(tj)).Lc = T_Lc(2);                                 % critical L for this track
    CI = [CI mean(L(:,2))];                                         % average confinement index per track
    MeanSqD = [MeanSqD; trajectory(a(tj)).msd(1:minimum)];    % msd first [minimum] time lags in um2
    clear Res 
    caption = sprintf('tj is %d from %d', tj, length(a));
  % Print to command window.
    fprintf('%s\n', caption);
%     L_all = [L_all L(:,2)];
end

% add running D for the remaining tracks (can later be used in synaptic
% script)
if strcmp(selected_tracks,'mobile') || strcmp(selected_tracks,'immobile')
    for tj = 1:length(b)
        xy = trajectory(b(tj)).xy;              % in um
        fr = trajectory(b(tj)).frame;           % framenumbers in track
        A = [xy fr];
        Dr = CalculateDr(A,dt,length(xy),npointsMSDset); % calculate running D in um2/s
        trajectory(b(tj)).Dr = Dr; 
    end
end

confinement_zones = [cz_index cz_x cz_y cz_rad Di cz_dwell cz_trans];               % used in dist_confinement_to_psd analysis
confinement_zones_transient = confinement_zones(cz_trans,:);
confined_index = find(confinement_data(:,4));
confined_data = confinement_data(confined_index,:);     % only trajectories with confinement zone(s)
fraction_confined = size(confined_data,1)/size(confinement_data,1);

if plotscatter
    % plot max difcoef vs confinement radius -> bigger difcoef, bigger radius can be found
    figure
    scatter(confined_data(:,13),confined_data(:,7))
    axis tight
    xlabel('Max diffusion coefficient [um2/s]');
    ylabel('Confinement radius [nm]');

    % plot track length vs confinement radius
    figure
    scatter(confined_data(:,2),confined_data(:,7))
    axis tight
    xlabel('Track length');
    ylabel('Confinement radius [nm]');

    % plot track length vs Dmax
    figure
    scatter(confined_data(:,2),confined_data(:,13))
    axis tight
    xlabel('Track length');
    ylabel('Max diffusion coefficient [um2/s]');
end

% SAVING
if savefile
    save('trajectories.mat','trajectory', '-append')
    save confinement.mat CI MeanSqD minimum confinement_data confined_data ...
        Sm Lcm Tc alfa Lmax Lcmax renderpx selected_tracks ...
        confinement_zones L_all fraction_confined border borderpx np ...
        npointsMSDset Smin confinement_zones_transient ...
        imm_conf_map radius_imm

    file_name = sprintf('confinement_data_%s.txt',selected_tracks);
    save(file_name, 'confinement_data', '-tabs', '-ascii')
    file_name = sprintf('confined_data_%s.txt',selected_tracks);
    save(file_name, 'confined_data', '-tabs', '-ascii')
end

% PLOT MAP OF CONFINEMENT ZONES OVERLAYED WITH TRACKS
% load trajectories
% load confinement
% extract_confined_param  % load parameters

index_tracklength = find([trajectory.length] >= minimum);

% which trajectories to analyze?
if strcmp(selected_tracks,'mobile')
    a = mobile;
elseif strcmp(selected_tracks,'immobile')
    a = immobile;
elseif strcmp(selected_tracks,'all')
    a = index_tracklength;
else
    error('Error: Define which tracks to analyze.')
end

xy = [];
for n=1:length(trajectory(a))
    xy = [xy ; trajectory(a(n)).xy];       %to get 2 columns with xy coordinates of all tracks below each other
end

xy = xy.*1000;  % from um to nm
x = xy(:,1); y = xy(:,2);

xx = ceil(x ./renderpx)+borderpx;               % +border to be able to plot 2d gaussian at the borders
yy = ceil(y ./renderpx)+borderpx;               % +border to be able to plot 2d gaussian at the borders
im = zeros(max(yy)+borderpx, max(xx)+borderpx); % +border to be able to plot 2d gaussian at the borders

% plot heatmap for tracks on coordinates
% for n = 1:length(xx)
%     im(yy(n), xx(n)) = im(yy(n), xx(n)) + 1;  % matrix with value 1 or higher when there is a trajectory on that coordinate
% end
% RI = imref2d(size(im));
% RI.XWorldLimits = [(0-(borderpx*renderpx)) (size(im,2)*renderpx-(borderpx*renderpx))];
% RI.YWorldLimits = [(0-(borderpx*renderpx)) (size(im,1)*renderpx-(borderpx*renderpx))];
% figure; imshow(im, RI, 'InitialMagnification','fit', 'Colormap', jet, 'DisplayRange', [0 5]);
% axis on; xlabel('x (nm)'); ylabel('y (nm)');

S = size(im);
confinement_map = zeros(S);
%a = find([trajectory.length]);

for tj = 1:length(a)
    Res = trajectory(a(tj)).R;
    cregions = size(Res);
    if cregions(1) >0    
        for c = 1:cregions(1)
            xc = ceil(Res(c,2)./renderpx)+borderpx;     % x
            yc = ceil(Res(c,3)./renderpx)+borderpx;     % y
            rc = ceil(Res(c,1)./renderpx);              % radius
            if rc > 0 && Res(c,2) > 0 && Res(c,3) > 0  % only confinement zones with positive center will be plotted
                % plot each zone as a 2d gauss with FWHM = domain radius (rc)
                % n.b. sigma = FWHM/2.355
                md = 2 * floor((rc*3)/2) + 1;           % create matrix with size rc*3 and round to nearest odd number
                m = zeros(md, md);
                cp = ceil(md/2);                        % center point
                mat = gauss2d(m, ceil(rc)/2.355, [cp cp]);
                confinement_map(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1)) = confinement_map(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1))+mat;
            end
        end
    end
end

% plot centers from immobile tracks as a heatmap in addition to confinement
% zones
if imm_conf_map
    a = cntr_immobile;
    for tj = 1:length(a)
        xc = ceil(a(tj,2).*1000./renderpx)+borderpx;    % x center from um to nm
        yc = ceil(a(tj,3).*1000./renderpx)+borderpx;    % y center from um to nm
        rc = ceil(radius_imm./renderpx);                      % radius
        if xc < max(xx) && yc < max(yy)
            if rc > 0 && a(tj,2) > 0 && a(tj,3) > 0  % only with positive center will be plotted
                % plot each zone as a 2d gauss with FWHM = domain radius (rc)
                % n.b. sigma = FWHM/2.355
                md = 2 * floor((rc*3)/2) + 1; % round to nearest odd number
                m = zeros(md, md);
                cp = ceil(md/2); % center point
                mat = gauss2d(m, ceil(rc)/2.355, [cp cp]);
                confinement_map(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1)) = confinement_map(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1))+mat;
            end
        end
    end
end

% plot heatmap of confinement zones
h = figure;
RI = imref2d(size(confinement_map));
RI.XWorldLimits = [(0-border) (size(confinement_map,2)*renderpx-border)];
RI.YWorldLimits = [(0-border) (size(confinement_map,1)*renderpx-border)];
load Cmap_confinement 
% imshow(confinement_map, 'Colormap', mycmap, 'DisplayRange', [0 2.5]); hold on;
imshow(confinement_map, RI, 'InitialMagnification','fit', 'Colormap', big_cmap, 'DisplayRange', [0 2.5]); hold on;
axis on; xlabel('x (nm)'); ylabel('y (nm)'); box off; colorbar; axis equal;

if savefile 
    file_name = sprintf('Confinement_map_%s',selected_tracks);
    saveas(gcf,file_name);
end

% plot tracks and confinement circles over heatmap
a = confined_data(:,1);
CM = jet(length(a));
for tj = 1:length(a)
    xy = trajectory(a(tj)).xy.*1000;  % from um to nm
    plot(xy(:,1), xy(:,2),'color', CM(tj,:));
%     text(trajectory(a(tj)).xy(1,1).*1000+(10),trajectory(a(tj)).xy(1,2).*1000, ...
%         num2str(a(tj)),'Color','k','FontSize',10);
    Res = trajectory(a(tj)).R;
    cregions = size(Res);
    if cregions(1) >0
        if ismember(a(tj),mobile) == true
            col = 'red';
        else
            col = 'black';
        end
        for c = 1:cregions(1)
            xc = Res(c,2);      % x
            yc = Res(c,3);      % y
            rc = Res(c,1);      % radius
            circle_color(xc, yc, rc, col);
        end
    end
end
axis equal;
% PlotBoundaries(filename_wf,pixel)
if savefile
    file_name = sprintf('Confinement_map_%s_circles',selected_tracks);
    saveas(gcf,file_name);
end
disp('done... data stored');

%% plot confinement map transient confinement and immobile tracks in separate figures
% PLOT HEATMAP OF CONFINEMENT ZONES
load trajectories
load confinement
extract_confined_param  % load parameters

a = mobile;

xy = [];
for n=1:length(trajectory)
    xy = [xy ; trajectory(n).xy];       %to get 2 columns with xy coordinates of all tracks below each other
end

xy = xy.*1000;  % from um to nm
x = xy(:,1); y = xy(:,2);

xx = ceil(x ./renderpx)+borderpx;               % +border to be able to plot 2d gaussian at the borders
yy = ceil(y ./renderpx)+borderpx;               % +border to be able to plot 2d gaussian at the borders
im = zeros(max(yy)+borderpx, max(xx)+borderpx); % +border to be able to plot 2d gaussian at the borders

S = size(im);
confinement_map1 = zeros(S);

% plot transient confinement zones
for tj = 1:length(a)
    Res = trajectory(a(tj)).R;
    cregions = size(Res);
    if cregions(1) >0    
        for c = 1:cregions(1)
            xc = ceil(Res(c,2)./renderpx)+borderpx;     % x
            yc = ceil(Res(c,3)./renderpx)+borderpx;     % y
            rc = ceil(Res(c,1)./renderpx);              % radius
            if rc > 0 && Res(c,2) > 0 && Res(c,3) > 0  % only confinement zones with positive center will be plotted
                % plot each zone as a 2d gauss with FWHM = domain radius (rc)
                % n.b. sigma = FWHM/2.355
                md = 2 * floor((rc*3)/2) + 1;           % create matrix with size rc*3 and round to nearest odd number
                m = zeros(md, md);
                cp = ceil(md/2);                        % center point
                mat = gauss2d(m, ceil(rc)/2.355, [cp cp]);
                confinement_map1(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1)) = confinement_map1(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1))+mat;
            end
        end
    end
end

% plot heatmap of confinement zones
h1 = figure;
RI = imref2d(size(confinement_map1));
RI.XWorldLimits = [(0-border) (size(confinement_map1,2)*renderpx-border)];
RI.YWorldLimits = [(0-border) (size(confinement_map1,1)*renderpx-border)];
load Cmap_confinement 
% imshow(confinement_map, 'Colormap', mycmap, 'DisplayRange', [0 2.5]); hold on;
imshow(confinement_map1, RI, 'InitialMagnification','fit', 'Colormap', big_cmap, 'DisplayRange', [0 2.5]); hold on;
axis on; xlabel('x (nm)'); ylabel('y (nm)'); box off; colorbar; axis equal;
a1 = gca;
if savefile 
    file_name = sprintf('Confinement_map_transient');
    saveas(gcf,file_name);
end

% a = immobile;
% 
% xy = [];
% for n=1:length(trajectory(a))
%     xy = [xy ; trajectory(a(n)).xy];       %to get 2 columns with xy coordinates of all tracks below each other
% end
% 
% x = xy(:,1); y = xy(:,2);
% 
% xx = ceil(x ./renderpx)+borderpx;               % +border to be able to plot 2d gaussian at the borders
% yy = ceil(y ./renderpx)+borderpx;               % +border to be able to plot 2d gaussian at the borders
% im = zeros(max(yy)+borderpx, max(xx)+borderpx); % +border to be able to plot 2d gaussian at the borders
% 
% S = size(im);
confinement_map2 = zeros(S);

% plot centers from immobile tracks as a heatmap
b = cntr_immobile;
for tj = 1:length(b)
    xc = ceil(b(tj,2).*1000./renderpx)+borderpx;    % x
    yc = ceil(b(tj,3).*1000./renderpx)+borderpx;    % y
    rc = ceil(radius_imm./renderpx);                      % radius
    if xc < max(xx) && yc < max(yy)
        if rc > 0 && b(tj,2) > 0 && b(tj,3) > 0  % only confinement zones with positive center will be plotted
            % plot each zone as a 2d gauss with FWHM = domain radius (rc)
            % n.b. sigma = FWHM/2.355
            md = 2 * floor((rc*3)/2) + 1; % round to nearest odd number
            m = zeros(md, md);
            cp = ceil(md/2); % center point
            mat = gauss2d(m, ceil(rc)/2.355, [cp cp]);
            confinement_map2(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1)) = confinement_map2(round(yc)-(cp-1):round(yc)+(cp-1), round(xc)-(cp-1):round(xc)+(cp-1))+mat;
        end
    end
end

% plot heatmap of confinement zones
h2 = figure;
RI = imref2d(size(confinement_map2));
RI.XWorldLimits = [(0-border) (size(confinement_map2,2)*renderpx-border)];
RI.YWorldLimits = [(0-border) (size(confinement_map2,1)*renderpx-border)];
load Cmap_confinement 
% imshow(confinement_map, 'Colormap', mycmap, 'DisplayRange', [0 2.5]); hold on;
imshow(confinement_map2, RI, 'InitialMagnification','fit', 'Colormap', big_cmap, 'DisplayRange', [0 2.5]); hold on;
axis on; xlabel('x (nm)'); ylabel('y (nm)'); box off; colorbar; axis equal;
a2 = gca;
if savefile 
    file_name = sprintf('Confinement_map_immobile');
    saveas(gcf,file_name);
end

linkaxes([a1 a2]);

disp('done... data stored');