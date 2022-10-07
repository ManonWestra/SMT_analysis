%% Single-molecule tracking analysis
% read localization csv file from ONI nanoimager, estimate diffusion
% coefficient and discriminate between mobile and immobile trajectories
%   (before exporting csv file in ONI software, check all boxes under 
%    Settings - User Settings - Export options)
% INPUT
% *.csv                 localization table (check the order of the required
%                       columns, might be different in other versions of the software)
% max_roi.tif           maximum projection widefield/tirf 
%                       can be used as mask
% max_roi_homer.tif     maximum projection of PSD marker
% OUTPUT
% create a data structure 'trajectory':
% for each trajectory the following params are stored in this structure:
% x-y       x and y coordinates [um]
% length    track length [# frames]
% frame     frame #
% MSD       MSD [um2]
% time      time [sec]
% Dinst     instantanous D [um2/s]
% mobile    whether particle is mobile based on ratio between radius of
%           gyration and mean stepsize
% immobile  whether particle is immobile based on ratio between radius of
%           gyration and mean stepsize
%
% Manon Westra, 2022
    
clear; clc;

convertResults2Trjectory_param  % load parameters

% INPUT FILES
filelist = dir(filename_csv);
file = filelist.name;
data = csvread(file,1,0);

if ROI_image == 0
    % create or load max projection image widefield
    if maxroi == 0
        imagelist = dir(filename_wf_raw);
        fileinfo = imfinfo(imagelist.name);
        imsizerow = fileinfo(1).Height; 
        imsizecol = fileinfo(1).Width;
        imsize = [imsizerow imsizecol];
        frames = size(fileinfo,1);
        widefield = MaxProjection(imagelist.name, imsize, frames);
        imwrite(widefield, filename_wf);
    elseif maxroi == 1
        imagelist = dir(filename_wf);
        widefield = imread(imagelist.name);
    end

    % create or load max projection image PSD
    if maxroi == 0
        imagelist = dir(filename_psd_raw);
        fileinfo = imfinfo(imagelist.name);
        imsizerow = fileinfo(1).Height; 
        imsizecol = fileinfo(1).Width;
        imsize = [imsizerow imsizecol];
        frames = size(fileinfo,1);
        homer = MaxProjection(imagelist.name, imsize, frames);
        imwrite(homer, filename_psd);
    elseif maxroi == 1
        imagelist = dir(filename_psd);
        homer = imread(imagelist.name);
    end
end

% ANALYSIS

load param
param.mem = disappear_steps;
disp(['loading ' file]);

% select channel for analysis
if channel == 1
    C = find(data(:,1));
elseif channel == 0
    C = find(~data(:,1));
else
    error('Error: Define which channel to track.')    
end

data = data(C,:);

% select localizations with localization precision <= loc_precision
err = (data(:,x_prec_nm)+ data(:,y_prec_nm))./2;    % X precision (nm) Y precision (nm)
S = find(err <= loc_precision);
data = [data(S,x_px) data(S,y_px) data(S, frame_num)]; % [x(px) y(px) frame]

% select localizations in rectangular roi [comment out if you analyze whole FOV]
%data = get_roi_coords(data);
%saveas(gcf,'rectangle.png');

% plot density map of the localizations
plot_density_map(data(:,1), data(:,2), pixel, rendered_pix, 1);
if savefile
    export_fig Density_map.png
end

if ROI_image == 1   % draw ROI on rendered image
    if load_prev_ROI
        % instead of drawing, load a previously drawn roi
        load roimask
    else
        roi_outline = drawfreehand(gca, 'MultiClick', true);
        if savefile
            saveas(gcf,'roi_selected');
        end
        roi_outline.Position(:,1) = roi_outline.Position(:,1).*(rendered_pix/pixel);
        roi_outline.Position(:,2) = roi_outline.Position(:,2).*(rendered_pix/pixel);
    end
    xroi = roi_outline.Position(:,1);
    yroi = roi_outline.Position(:,2);
    inoutroi = inpolygon (data(:,1), data(:,2), xroi, yroi);
elseif ROI_image == 0  % draw ROI on widefield image
    if load_prev_ROI
        % instead of drawing, load a previously drawn roi
        load roimask
    else
%         widefield = imread(filename_wf);
        [M, N] = size(widefield);
        % image to world coordinates (first pixel 0,0 instead of 1,1)
        RI = imref2d(size(widefield));
        RI.XWorldLimits = [(-0.5) (size(widefield,2)-(0.5))];
        RI.YWorldLimits = [(-0.5) (size(widefield,1)-(0.5))];
        % select localizations inside ROI from widefield image for analysis
        figure; imshow(widefield,RI,[lower_limit_image upper_limit_image]);
        roi_outline = drawfreehand(gca, 'MultiClick', true);
        if savefile
            saveas(gcf,'roi_selected');
        end
    end
    xroi = roi_outline.Position(:,1);
    yroi = roi_outline.Position(:,2);
    inoutroi = inpolygon (data(:,1), data(:,2), xroi, yroi);
end
% close
if savefile
    save('roimask.mat', 'roi_outline');
end

data = data(inoutroi,:);

% tracking
data = sortrows(data,3);                        % sort on framenumber
tracks = track(data, tracking_radius, param);   % results in [x y frame unique_ID_number]
clc;

% convert data to nm 
tracks(:, 1:2) = tracks(:, 1:2) .* pixel;

% Generate a list of lengths of each track
[t_lengths,track_ids] = ...
hist(tracks(:, 4),unique(tracks(:, 4)));
% find all tracks in the range of track_min_max
tracks_long = track_ids((t_lengths >= track_min_max(1)) & (t_lengths <= track_min_max(2)));

% Plot the histogram of the lengths with bins from 1 to the number
% of frames in each segment
figure
hist(t_lengths(tracks_long),1:max(t_lengths(tracks_long)))
xlabel('Track length'); ylabel('Frequency');
if savefile
    saveas(gcf,'Length_histogram.png');
end
close

ntracks = length(tracks_long);         % # tracks
disp(['analyzing ' num2str(ntracks) ' tracks...']);

% analyze individual tracks and combine data in structure 'trajectory'
dif_coef = [];
trajectory_length = zeros(ntracks,1);

for tj = 1:length(tracks_long)
    ff = tracks(:, 4) == tracks_long(tj);    % find localizations with same track ID as tj
    ff = find(ff);
    xy = tracks(ff, 1:2)./1000;     % x y coordinates in um 
    A = [xy tracks(ff,3)];  % x y frame
    % store data
    trajectory_length(tj) = length(tracks(ff, 1:2));    % number of localizations with same track ID
    trajectory(tj).xy = xy;  % x y coordinates in um
    trajectory(tj).length = length(tracks(ff, 1:2));
    trajectory(tj).frame = tracks(ff, 3); 
    trajectory(tj).msd = MSD_new(A,dt);    % [MSD stdev n timelag] (works if frames are missing) um2
    trajectory(tj).time = tracks(ff,3)*dt;
    % fit MSD to estimate Dinst
    if length(xy) >= np
        D = Dinst(trajectory(tj).msd(:,4),trajectory(tj).msd(:,1),np);       % estimate instantaneous D in um2/s
        trajectory(tj).Dinst = D;
        dif_coef = [dif_coef; D];
    end
        
  caption = sprintf('tj is %d from %d', tj, ntracks);
  % Print to command window.
  fprintf('%s\n', caption);
end

% plot all tracks in random colors
figure('Color', 'white'); hold on;
for tj = 1:length(trajectory) 
    plot(trajectory(tj).xy(:,1), -trajectory(tj).xy(:,2)) % minus before y values to get same orientation as image
end
axis equal; hold off;
xlabel('x (um)'); ylabel('y (um)');
if savefile
    saveas(gcf,'Trajectory_map_all');
    export_fig Trajectory_map_all.pdf
end

% discriminate between mobile/immobile particles based on gyration radius and stepsize
NR = [];
Dm = [];
Di = [];
mobile = [];     % can be used as index for mobile trajectories
immobile = [];   % can be used as index for immobile trajectories
cntr_immobile = [];  % store centers of immobile tracks
figure('Color', 'white');hold on;

for tj = 1:length(trajectory)
    xy = trajectory(tj).xy;
    N = length(xy);
    mx = mean(xy(:, 1));
    my = mean(xy(:, 2));
    Rg = sqrt(1/N * sum((xy(:,1) - mx).^2 + (xy(:,2)-my).^2));
    deltaCoords = xy(2:end, :) - xy(1:end-1, :);
    squaredDisplacement = sum(deltaCoords.^2, 2); %# dx^2+dy^2
    mean_step_size = sqrt(mean(squaredDisplacement));
    norm_ratio = (sqrt(pi/2)*Rg)/mean_step_size;
           
    NR = [NR; norm_ratio];
    D = trajectory(tj).Dinst;
    
    if norm_ratio >= mob_ratio
        col = 'red';
        Dm = [Dm; D];   % diffusion coefficient of mobile trajectories
        mobile = [mobile; tj];
        trajectory(tj).mobile = true;
    else
        col = 'black';
        Di = [Di; D];   % diffusion coefficient of immobile trajectories
        immobile = [immobile; tj];
        trajectory(tj).immobile = true;
        cntr_immobile = [cntr_immobile;tj mx my];
    end
    plot(xy(:,1), -xy(:,2), 'color', col);
end
axis equal; axis off; hold off;
xlabel('x (um)'); ylabel('y (um)');
if savefile
    export_fig immobile_mobile_trjectories.png
end 
close
% plot diffusion coefficient mobile/immobile log D [um^2/s]
figure('Color', 'white'); hold on;
x = linspace(-5, 2,50);
y1 = hist(log10(Dm), x); y1=(y1)./sum(y1); 
y2 = hist(log10(Di), x); y2=(y2)./sum(y2);
plot(x, y1, 'red'); plot(x, y2,'black');
title('Instantaneous diffusion coefficient immobile/mobile'); set(gca,'XLim',[-5 2]);xlabel('log D [um^2/s]');
hold off;
if savefile
    saveas(gcf,'Dinst_mob_immob.png');
end 
close
% fraction of mobile particles of total longer than minSteps
fraction_mobile = size(Dm,1)/size(trajectory,2);

% trajectory struct without immobile particles
%trajectory = trajectory(mobile);
%trajectory_length = trajectory_length(mobile);

% plot tracks in random colors
figure('Color', 'white'); hold on;
for tj = 1:length(trajectory)   
    plot(trajectory(tj).xy(:,1), -trajectory(tj).xy(:,2))
end
axis equal; axis off; hold off;
xlabel('x (um)'); ylabel('y (um)');
if savefile
    saveas(gcf,'Trajectory_map');
    export_fig Trajectory_map.pdf
end
close
% diffusion coefficient from D [um^2/s] to log D [um^2/s]
Deff = log10(dif_coef);
Deff_imm = log10(Di);
Deff_mob = log10(Dm);

if savefile
    save trajectories.mat trajectory trajectory_length dif_coef ...
        Deff pixel dt tracking_radius minSteps minPlot ...
        mobile mob_ratio immobile fraction_mobile disappear_steps channel ...
        Deff_imm Deff_mob NR loc_precision track_min_max rendered_pix np ...
        cntr_immobile
end
disp('done... data stored');

% plot diffusion coefficient of all trajectories
figure('Color', 'white'); hold on;
x = linspace(-5, 2,100);
y = hist(Deff, x); y=y./max(y);
bar(x,y)
title('Instantaneous diffusion coefficient'); set(gca,'XLim',[-5 2]);xlabel('log D [um^2/s]');set(gca,'YLim',[0 1]);
hold off;
if savefile
    saveas(gcf,'Dinst.png');
end
% close

if ROI_image == 0  % plot only if there are widefield images
    % image to world coordinates (first pixel 0,0 instead of 1,1)
    RI = imref2d(size(widefield));
    RI.XWorldLimits = [(-(0.5*pixel/1000)) (size(widefield,2)*pixel/1000-(0.5*pixel/1000))];
    RI.YWorldLimits = [(-(0.5*pixel/1000)) (size(widefield,1)*pixel/1000-(0.5*pixel/1000))];

    % plot all trajectories longer than [minPlot] on widefield image
    figure;
    imshow(widefield,RI,[lower_limit_image upper_limit_image]);  %, [300 700]
    hold on;
    CM = jet(length(trajectory));
    for tj = 1:length(trajectory)
        plot(trajectory(tj).xy(:,1), trajectory(tj).xy(:,2),'color', CM(tj,:))
    end
    hold off
    xlabel('x (um)'); ylabel('y (um)');
    a1 = gca;
    if savefile
        export_fig Tracks_on_widefield.png
    end

    % plot all trajectories longer than [minPlot] on homer
    % homer = imread(filename_psd);
    figure;
    imshow(homer,RI, [low_homer up_homer]);
    hold on;
    for tj = 1:length(trajectory)
        plot(trajectory(tj).xy(:,1), trajectory(tj).xy(:,2),'color', CM(tj,:))
    end
    hold off
    xlabel('x (um)'); ylabel('y (um)');
    a2 = gca;
    if savefile
        export_fig Tracks_on_homer.png
    end

    linkaxes([a1 a2]); % to zoom similarly on both images
end
