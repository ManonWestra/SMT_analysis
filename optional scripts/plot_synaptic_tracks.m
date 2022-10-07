%% Plot synaptic tracks over synaptic marker
load trajectories.mat
load PSDmask.mat
plot_synaptic_tracks_param  % load parameters

%Load max projection image
widefield = imread(filename_psd);
% image to world coordinates (first pixel 0,0 instead of 1,1)
RI = imref2d(size(widefield));
RI.XWorldLimits = [(-(0.5*pixel)) (size(widefield,2)*pixel-(0.5*pixel))];
RI.YWorldLimits = [(-(0.5*pixel)) (size(widefield,1)*pixel-(0.5*pixel))];

figure;
imshow(widefield, RI, [1 1000]);
hold on;

a =[];
for tj = 1:length(trajectory)
    if (trajectory(tj).Dsynaptic > 0)
        a = [a; tj];
    end
end

%figure('Color', 'white'); hold on;

tracks = trajectory(a);

for tj = 1:length(tracks)   
    plot(tracks(tj).xy(:,1)*1000, tracks(tj).xy(:,2)*1000)
end
axis equal; axis off; hold off;

%% Plot extrasynaptic tracks over synaptic marker
load trajectories.mat
load PSDmask.mat
plot_synaptic_tracks_param  % load parameters

%Load max projection image
widefield = imread(filename_psd);
% image to world coordinates (first pixel 0,0 instead of 1,1)
RI = imref2d(size(widefield));
RI.XWorldLimits = [(-(0.5*pixel)) (size(widefield,2)*pixel-(0.5*pixel))];
RI.YWorldLimits = [(-(0.5*pixel)) (size(widefield,1)*pixel-(0.5*pixel))];

figure;
imshow(widefield, RI, [1 1000]);
hold on;

a =[];
for tj = 1:length(trajectory)
    if (trajectory(tj).Dextrasynaptic > 0) 
        a = [a; tj];
    end
end

%figure('Color', 'white'); hold on;

tracks = trajectory(a);

for tj = 1:length(tracks)   
    plot(tracks(tj).xy(:,1)*1000, tracks(tj).xy(:,2)*1000)
end
axis equal; axis off; hold on;

%% Plot all tracks over synaptic marker in different color according to extra or synaptic

% Load max projection image
% widefield = imread('max_roi_homer.tif');
% % image to world coordinates (first pixel 0,0 instead of 1,1)
% RI = imref2d(size(widefield));
% RI.XWorldLimits = [(-(0.5*pixel)) (size(widefield,2)*pixel-(0.5*pixel))];
% RI.YWorldLimits = [(-(0.5*pixel)) (size(widefield,1)*pixel-(0.5*pixel))];
% 
% figure;
% imshow(widefield, RI, [1 1000]);
openfig(filename_psdmask);
hold on;

tracks = trajectory;

for tj = 1:length(tracks)
    if tracks(tj).Dsynaptic > 0
        plot(tracks(tj).xy(:,1)*1000, tracks(tj).xy(:,2)*1000, 'color', 'red')
    else
        plot(tracks(tj).xy(:,1)*1000, tracks(tj).xy(:,2)*1000, 'color', 'black')
    end
        
end
axis equal; hold off;     

saveas(gcf, 'PSD_with_synaptictracks');
export_fig Trajectory_map_synaptic.pdf

%% Plot all tracks in different color according to extra or synaptic

figure('color', 'white');hold on;

tracks = trajectory;

for tj = 1:length(tracks)
    if tracks(tj).Dsynaptic > 0
        plot(tracks(tj).xy(:,1), -tracks(tj).xy(:,2), 'color', 'red')
    else
        plot(tracks(tj).xy(:,1), -tracks(tj).xy(:,2), 'color', 'black')
    end
        
end
axis equal; axis off; hold off;
