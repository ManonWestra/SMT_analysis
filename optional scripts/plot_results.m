%% collection of different sections to plot results
% all based on results from other single-molecule tracking analysis scripts
% from:
%   convertResults2Trjectory
%   extract_confined
%   extract_extra_transient_synaptic_tracks
%   make_psd_mask
%
% Manon Westra, 2022

%% plot all trajectories
load trajectories

figure;hold on;
for tj = 1:length(trajectory)
    plot(trajectory(tj).xy(:,1), -trajectory(tj).xy(:,2)) % - before y values to get same orientation as image
end
hold off
axis equal
xlabel('x (um)'); ylabel('y (um)');

%% plot all trajectories longer than [minPlot] on PSD mask in nm
load trajectories
convertResults2Trjectory_param

openfig(filename_psdmask); 
hold on;
for tj = 1:length(trajectory)
    plot(trajectory(tj).xy(:,1).*1000, trajectory(tj).xy(:,2).*1000)
end
hold off
xlabel('x (nm)'); ylabel('y (nm)');
saveas(gcf,'Tracks_on_PSDmask');

%% plot all trajectories longer than [minPlot] on PSD mask in nm with track
% numbers next to it
load trajectories
convertResults2Trjectory_param

openfig(filename_psdmask); 
hold on;
for tj = 1:length(trajectory)
    plot(trajectory(tj).xy(:,1).*1000, trajectory(tj).xy(:,2).*1000)
    text(trajectory(tj).xy(1,1).*1000+(10),trajectory(tj).xy(1,2).*1000, num2str(tj),'Color','r',...
            'FontSize',10);
end
hold off
xlabel('x (nm)'); ylabel('y (nm)');

%% plot confinement zones over widefield image
load trajectories.mat
load confinement.mat
extract_confined_param  % load parameters

figure; hold on;
widefield = imread(filename_wf);

% image to world coordinates (first pixel 0,0 instead of 1,1)
RI = imref2d(size(widefield));
RI.XWorldLimits = [(-(0.5*pixel)) (size(widefield,2)*pixel-(0.5*pixel))];
RI.YWorldLimits = [(-(0.5*pixel)) (size(widefield,1)*pixel-(0.5*pixel))];

imshow(widefield, RI, [1 1000]);
hold on;
a = confined_data(:,1);
for tj = 1:length(a)
    xy = trajectory(a(tj)).xy;
    plot(xy(:,1),xy(:,2));
    Res = trajectory(a(tj)).R;
    cregions = size(Res);
    if cregions(1) >0
        if ismember(a(tj),mobile) == true
            col = 'red';
        else
            col = 'black';
        end
        for c = 1:cregions(1)
            xc = (Res(c,2));   % x
            yc = (Res(c,3));   % y
            rc = (Res(c,1));      % radius
            circle_color(xc, yc, rc, col);
        end
    end
end
axis equal;

%% plot filled PSD boundaries on top of confinement map
load PSDmask.mat
load trajectories.mat
extract_confined_param  % load parameters

openfig(filename_confmap);
for q = 1:length(B) 
    xc = B{q}(:,1);
    yc = B{q}(:,2);
    g = fill(yc, xc,'black'); set(g,'facealpha',.2);set(g,'LineStyle', 'none')   %transparent black PSDs
end

%% plot confinement zones over PSD mask
%figure; hold on;
load trajectories
load confinement
extract_confined_param  % load parameters
openfig(filename_psdmask); hold on; 
a =confined_data(:,1);    %synaptic_tracks
for tj = 1:length(a)
    xy = trajectory(a(tj)).xy;
    plot(xy(:,1),xy(:,2));
    Res = trajectory(a(tj)).R;
    cregions = size(Res);
    if cregions(1) >0
        if ismember(a(tj),mobile) == true
            col = 'red';
        else
            col = 'black';
        end
        for c = 1:cregions(1)
            xc = (Res(c,2));   % x
            yc = (Res(c,3));   % y
            rc = (Res(c,1));      % radius
            circle_color(xc, yc, rc, col);
        end
    end
end
axis equal;

%% plot all trajectories PSD mask in nm next to plot with synaptic/extrasynaptic track
openfig(filename_psdmask); 
hold on;
for tj = 1:length(trajectory)
    plot(trajectory(tj).xy(:,1).*1000, trajectory(tj).xy(:,2).*1000)
end
hAxes = findobj(gcf, 'type', 'axes');
hold off
xlabel('x (nm)'); ylabel('y (nm)');
saveas(gcf,'Tracks_on_PSDmask');

% Plot all tracks over synaptic marker in different color according to extra or synaptic
openfig(filename_psdmask);
hold on;

tracks = trajectory;
for tj = 1:length(tracks)
    if tracks(tj).Dsynaptic > 0
        plot(tracks(tj).xy(:,1).*1000, tracks(tj).xy(:,2).*1000, 'color', 'red')
    else
        plot(tracks(tj).xy(:,1).*1000, tracks(tj).xy(:,2).*1000, 'color', 'black')
    end
        
end
axis equal; 
hAxes2 = findobj(gcf, 'type', 'axes'); 
hold off;     
saveas(gcf, 'PSD_with_synaptictracks');

linkaxes([hAxes hAxes2])

%% MSD vs time plot & histogram Dinst
load trajectories.mat
load confinement.mat
% plot results
figure('Color', 'white');
% MSD vs time plot
% histogram trajectory length
subplot(1, 2, 1); 
x = 0:dt:dt*(minimum-1);
[r c] = size(MeanSqD);
A = mean(MeanSqD, 1); S = std(MeanSqD, 0, 1)/sqrt(r);
shadedErrorBar(x, A, S, 'lineProps', 'red');

title('average MSD vs time'); set(gca,'XLim',[0 minimum*dt],'YLim',[0 0.5]);xlabel('interval (s)');
ylabel('MSD (um^2)');

% histogram Dinst
%figure('Color', 'white', 'Position', [900 1250 1200 400]);
x = linspace(-5, 2,100);
%Deff = log10(dif_long_track./10^6);
% Deff = log10(Deff);
y = hist(Deff, x); y=y./max(y);
subplot(1, 2, 2); 

title('instantaneous diffusion coefficient'); set(gca,'XLim',[-5 2]);xlabel('log D (um^2/s)');set(gca,'YLim',[0 1]);
ylabel('relative frequency (fractions)');
hold on;
%
options = statset('MaxIter',1000);
GMM = fitgmdist(Deff, 2,'Options',options);
GMM.mu
GMM.ComponentProportion

yfit1 = pdf('normal', x, GMM.mu(1), sqrt(GMM.Sigma(1)),r)*GMM.ComponentProportion(1);
yfit2 = pdf('normal', x, GMM.mu(2), sqrt(GMM.Sigma(2)),r)*GMM.ComponentProportion(2);
bar(x, y)
%
plot(x, yfit1);plot(x, yfit2);

%% plot track, D vs time and CI vs time for individual trajectories
load trajectories.mat
load confinement.mat
% n = tracknumber;
n = 114;
CM = jet(length(trajectory(n).xy(:,1)));

figure('Color','white');
%
hold on
scatter(trajectory(n).xy(1,1)*1000, -trajectory(n).xy(1,2)*1000, 'green','filled');
scatter(trajectory(n).xy(end,1)*1000, -trajectory(n).xy(end,2)*1000, 'red','filled');

for t = 1:length(trajectory(n).xy(:,1))-1
%     if trajectory(n).L(t)>31.6
%         plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'green');
%     else
%         plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'black');
%     end
%         
      plot(trajectory(n).xy(t:t+1,1)*1000, -trajectory(n).xy(t:t+1,2)*1000,'color',CM(t,:));
end
axis equal
% plot confinement circles (optional, comment out if not)
Res = trajectory(n).R;
    cregions = size(Res);
    if cregions(1) >0
            col = 'black';
        for c = 1:cregions(1)
            xc = Res(c,2);      % x
            yc = -Res(c,3);      % y
            rc = Res(c,1);      % radius
            circle_color(xc, yc, rc, col);
        end
    end


figure('Color','white');

subplot(2,1,1)
hold on
for t = 1:length(trajectory(n).Dr)-1
    plot(trajectory(n).time(t:t+1), trajectory(n).Dr(t:t+1),'color',CM(t,:), 'LineWidth',1)
    %plot(time(t:t+1),trajectory(n).Dr(t:t+1),'color',CM(t,:))
end
title('Diffusion');xlabel('Time (s)');
% set(gca,'YLim',[0 0.3])
% set(gca,'XLim',[0 25])

subplot(2,1,2)
hold on
%scatter(trajectory(n).time(1:4), trajectory(n).L(1:4))
for t = 1:length(trajectory(n).L)-1
    plot(trajectory(n).time(t:t+1), trajectory(n).L(t:t+1),'color',CM(t,:), 'LineWidth',1)
    %plot(time(t:t+1),trajectory(n).L(t:t+1),'color',CM(t,:))
end
title('L');xlabel('Time (s)');
% set(gca,'XLim',[0 25])
% plot Lc threshold on top of L vs time plot
yline(trajectory(n).Lc,'LineStyle', '- -');
% thr_plot = 1:size(trajectory(n).time,2);
% thr_plot = thr_plot*trajectory(n).Lc./thr_plot;
% plot(trajectory(n).time,thr_plot, 'Color','black', 'LineStyle', '- -');

%% plot track and D vs time for individual trajectories marked for inout PSD
% n = tracknumber;
n = 114;
CM = jet(length(trajectory(n).xy(:,1)));

figure('Color','white');
%
hold on
scatter(trajectory(n).xy(1,1)*1000, -trajectory(n).xy(1,2)*1000, 'green','filled');
scatter(trajectory(n).xy(end,1)*1000, -trajectory(n).xy(end,2)*1000, 'red','filled');

for t = 1:length(trajectory(n).xy(:,1))-1
     if trajectory(n).inout(t)>0
         plot(trajectory(n).xy(t:t+1,1)*1000, -trajectory(n).xy(t:t+1,2)*1000,'red');
     else
         plot(trajectory(n).xy(t:t+1,1)*1000, -trajectory(n).xy(t:t+1,2)*1000,'black');
     end
end
axis equal

figure('Color','white');

subplot(2,1,1)
hold on
for t = 1:length(trajectory(n).Dr)-1
    if trajectory(n).inout(t)>0
        cc = 'red';
        lw = 2;
    else
        cc = 'black';
        lw = 1;
    end
    plot(trajectory(n).time(t:t+1), trajectory(n).Dr(t:t+1),'color',cc, 'LineWidth', lw)
end
title('diffusion');xlabel('time (s)');
% set(gca,'YLim',[0 1.2])

subplot(2,1,2)
hold on
%scatter(trajectory(n).time(1:4), trajectory(n).L(1:4))
for t = 1:length(trajectory(n).L)-1
    if trajectory(n).inout(t)>0
        cc = 'red';
        lw = 2;
    else
        cc = 'black';
        lw = 1;
    end
    plot(trajectory(n).time(t:t+1), trajectory(n).L(t:t+1),'color',cc, 'LineWidth', lw)
end
title('L');xlabel('time (s)');

%% Plot tracks and color-code for diffusion coefficient
figure('Color', 'white');hold on;

long_track = find(trajectory_length >= minSteps);
dif_long_track = [];
for tj = 1:length(long_track)
    Dinst = trajectory(long_track(tj)).Dinst./10^6;
    dif_long_track = [dif_long_track; Dinst];
end

col = jet(length(dif_long_track));
qmin = min(dif_long_track); 
qmax = max(dif_long_track);

for tj = 1:length(long_track)
    xy = trajectory(long_track(tj)).xy;
    q = round((((trajectory(long_track(tj)).Dinst./10^6) - qmin)/(qmax-qmin)).*(length(dif_long_track)-1))+1;
    plot(xy(:,1)./pixel, -xy(:,2)./pixel, 'color', col(q,:));
end
axis equal; axis off; colormap Jet

%% plot track and color-code for Dinst
plot_results_param  % load parameters
n = tracknumber;

figure('Color','white');
col = jet;
%
hold on
scatter(trajectory(n).xy(1,1), trajectory(n).xy(1,2), 'green','filled');
scatter(trajectory(n).xy(end,1), trajectory(n).xy(end,2), 'red','filled');

qmin = min(trajectory(n).Dr); qmax = max(trajectory(n).Dr);


for t = 1:length(trajectory(n).xy(:,1))-1
    q = round((((trajectory(n).Dr(t)) - qmin)/qmax).*200)+1;
    plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'color',col(q,:))
    %plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'color',CM(t,:))
end
axis equal; colormap Jet


%% plot single track
n = tracknumber;

figure('Color', 'white'); hold on;
for t = 1:length(trajectory(n).xy(:,1))-1
    plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'Color', 'black')
end
axis equal

%% plot single track with thicker line part
n = tracknumber;

figure('Color', 'white'); hold on;
for t = 1:length(trajectory(n).xy(:,1))-1
    if t < 30
        plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'Color', 'black')
    elseif t >= 30 && t < 60
        plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'Color', 'black','LineWidth', 2)
    else
        plot(trajectory(n).xy(t:t+1,1), trajectory(n).xy(t:t+1,2),'Color', 'black')
    end
end
axis equal

%% plot confined tracks on PSD mask with confinement zones and tracknumbers
a = confined_data(:,1);
tracks = trajectory(a);
openfig(filename_psdmask); hold on;
% figure('Color', 'white'); hold on;
for tj = 1:length(tracks)
plot(tracks(tj).xy(:,1)*1000, tracks(tj).xy(:,2)*1000)
    text(tracks(tj).xy(1,1)*1000+(10),tracks(tj).xy(1,2)*1000, num2str(a(tj)),'Color','r',...
            'FontSize',10);
%plot(tracks(tj).xy(:,1), tracks(tj).xy(:,2))
    Res = tracks(tj).R;
    cregions = size(Res);
    if cregions(1) >0
        if ismember(a(tj),mobile) == true
            col = 'black';
        else
            col = 'red';
        end
        for c = 1:cregions(1)
            xc = (Res(c,2));   % x
            yc = (Res(c,3));   % y
            rc = (Res(c,1));      % radius
            circle_color(xc, yc, rc, col);
        end
    end
end
axis equal; hold off;

%% plot trajectories; color-coded for above confinement threshold

figure('Color', 'white'); hold on;

% select tracks
%a = find(trajectory_length >= minSteps);
a = mobile;

tracks = trajectory(a);

for tj = 1:length(tracks)
    L = tracks(tj).L.';
    for t = 1:length(tracks(tj).L)-1
        if L(t) > Lcm
            cc = 'red';
        else
            cc = 'black';
        end
        plot(tracks(tj).xy(t:t+1,1), -tracks(tj).xy(t:t+1,2),  'color', cc)
    end
end
axis equal;hold off;

%% plot trajectories; color-coded for inout PSD

figure('Color', 'white'); hold on;

% select tracks
%a = find (trajectory_length > 500);
% a = tracknumber;
%tracks = trajectory(a);
tracks = trajectory;

for tj = 1:length(tracks)
    inout = tracks(tj).inout.';
    for t = 1:length(tracks(tj).inout)-1
        if inout(t) >0
            cc = 'red';
        else
            cc = 'black';
        end
        plot(tracks(tj).xy(t:t+1,1), -tracks(tj).xy(t:t+1,2), 'color', cc)
    end
end
axis equal;hold off;

%% plot all trajectories with dotted widefield outline
load trajectories

figure('Color', 'white'); hold on;

% select tracks
%a = find(trajectory_length >= minSteps);
a = mobile;

tracks = trajectory(a);

for tj = 1:length(tracks)
    plot(tracks(tj).xy(:,1), tracks(tj).xy(:,2))
end
PlotBoundaries('max_roi.tif',pixel/1000)
hold off
axis equal

%% plot all trajectories with track numbers next to it
load trajectories

figure;hold on;
for tj = 1:length(trajectory)
    plot(trajectory(tj).xy(:,1), -trajectory(tj).xy(:,2))
    text(trajectory(tj).xy(1,1)+(0.2),-trajectory(tj).xy(1,2), num2str(tj),'Color','k',...
            'FontSize',10);
end
hold off
axis equal
xlabel('x (um)'); ylabel('y (um)');