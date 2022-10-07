%% calculate the distance of confinement zones to border PSD for all tracks
% multiple sections to calculate different distances
% Manon Westra, 2022

clear;clc;
load trajectories.mat
load confinement.mat
load PSDmask.mat

cz_x = confinement_zones(:,2);  % x center confinement zone (nm)
cz_y = confinement_zones(:,3);  % y center confinement zone (nm)

confinement_distance = zeros(length(cz_x),2);
confinement_distance(:,1) = confinement_zones(:,1);  % trajectory number

% loop over all PSDs and check for the shortest distance
for k = 1:length(B)
    xv = B{k}(:,2); yv = B{k}(:,1);
    d = p_poly_dist(cz_x,cz_y,xv,yv);
    for cz = 1:length(cz_x)
        if d(cz) < confinement_distance(cz,2) || confinement_distance(cz,2) == 0
            confinement_distance(cz,2) = d(cz);
        end
    end
end

figure('color', 'white');
hist(confinement_distance(:,2),100);

% plot histogram confinement distance within 1000nm from border PSD
figure('color', 'white');
confinement_distance1000_ind = confinement_distance(:,2) <= 1000;
confinement_distance1000 = confinement_distance(confinement_distance1000_ind,:);
hist(confinement_distance1000(:,2),100);
percConfInPSD = sum(confinement_distance1000(:,2) <= 0)/length(confinement_distance1000)*100;
export_fig confinement_distance1000.png
save confinement_distance.mat confinement_distance confinement_distance1000 percConfInPSD

%% calculate distance immobile tracks centers to PSD border
clear;clc;
load trajectories.mat
load PSDmask.mat

cz_x = cntr_immobile(:,2).*1000;  % x center nm
cz_y = cntr_immobile(:,3).*1000;  % y center nm

imm_confinement_distance = zeros(length(cz_x),2);
imm_confinement_distance(:,1) = cntr_immobile(:,1);  % trajectory number

% loop over all PSDs and check for the shortest distance
for k = 1:length(B)
    xv = B{k}(:,2); yv = B{k}(:,1);
    d = p_poly_dist(cz_x,cz_y,xv,yv);
    for cz = 1:length(cz_x)
        if d(cz) < imm_confinement_distance(cz,2) || imm_confinement_distance(cz,2) == 0
            imm_confinement_distance(cz,2) = d(cz);
        end
    end
end

figure('color', 'white');
hist(imm_confinement_distance(:,2),100);

% plot histogram confinement distance within 1000nm from border PSD
figure('color', 'white');
imm_confinement_distance1000_ind = imm_confinement_distance(:,2) <= 1000;
imm_confinement_distance1000 = imm_confinement_distance(imm_confinement_distance1000_ind,:);
hist(imm_confinement_distance1000(:,2),100);
percImmInPSD = sum(imm_confinement_distance1000(:,2) <= 0)/length(imm_confinement_distance1000)*100;
export_fig imm_confinement_distance1000.png
save imm_confinement_distance.mat imm_confinement_distance imm_confinement_distance1000 percImmInPSD

%% distance confinement to PSD only for PSD associated tracks
clear;clc;
load trajectories.mat
load confinement.mat
load PSDmask.mat

synaptic = [trajectory.synaptic].';
synaptic_tracks = find(synaptic ~= 0);

ff = ismember(confinement_zones(:,1),synaptic_tracks);
confinement_zones_psd = confinement_zones(ff,:);
cz_x = confinement_zones_psd(:,2);
cz_y = confinement_zones_psd(:,3);

confinement_distance_psdtracks = zeros(length(cz_x),2);
confinement_distance_psdtracks(:,1) = confinement_zones_psd(:,1); % trajectory number

for k = 1:length(B)
    xv = B{k}(:,2); yv = B{k}(:,1);
    d = p_poly_dist(cz_x,cz_y,xv,yv);
    for cz = 1:length(cz_x)
        if d(cz) < confinement_distance_psdtracks(cz,2) || confinement_distance_psdtracks(cz,2) == 0
            confinement_distance_psdtracks(cz,2) = d(cz);
        end
    end
end

figure('color', 'white');
hist(confinement_distance_psdtracks(:,2),100);
export_fig confinement_distance_psdtracks.png

save('confinement_distance_psd.mat','confinement_distance_psdtracks')

%% works for single PSD analysis   from extract_extra_transient_synaptic_tracks script
% Results in tracks associated with only one specified PSD
dist_confinement_to_PSD_param   % load PSD number

synaptic = [trajectory.synaptic].';
synaptic_tracks = find(synaptic ~= 0);
tracks = trajectory(synaptic_tracks);

ff = ismember(confinement_zones(:,1),synaptic_tracks);
confinement_zones_psd = confinement_zones(ff,:);
cz_x = confinement_zones_psd(:,2);
cz_y = confinement_zones_psd(:,3);

load PSDmask.mat
xv = B{psd_number}(:,2);
yv = B{psd_number}(:,1);

% measure distances of confinement zones surrounding one specified PSD
distance_psd = [];

distance_psd = p_poly_dist(cz_x,cz_y,xv,yv);


