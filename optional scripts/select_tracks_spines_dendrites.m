%% select tracks in spines and/or dendrites
% based on ROIs made with 'draw_ROIs_on_image.m'
% output is the list of diffusion coefficients per group
% last section is plotting for comparison
% 
% Manon Westra, April 2021

clear;clc;
load trajectories.mat
load spinedendrite.mat

select_tracks_spines_dendrites_param  % load parameters

%% select tracks in dendrites
dendritetracks = NaN(length(trajectory),1);

% for each dendrite ROI find overlapping tracks
for k = 1:length(dendrites)
    xv = dendrites{k}(:,2); yv = dendrites{k}(:,1);
    for tj = 1:length(trajectory)
        if trajectory(tj).length >= MTL 
            inout = inpoly2(trajectory(tj).xy*1000, [xv yv]);
            dendrite = sum(inout)/trajectory(tj).length;   % fraction of track in dendrite
                       
            % overlapping with dendrite ROI
            if dendrite > 0 
                dendritetracks(tj) = dendrite;
            end
        end
    end        
end

% create structure with only dendrite associated tracks
dendrite_trajectory = trajectory(dendritetracks>0);

% list of diffusion coefficients in log(um2/s)
Ddendrite = (log10([dendrite_trajectory.Dinst]))';

% determine fraction of dendrite tracks
fraction_dendrite = length(Ddendrite)/length(trajectory);

if savefile
    save dendrite.mat MTL Ddendrite fraction_dendrite dendritetracks
    save Ddendrite.txt Ddendrite -tabs -ascii
end

%% select tracks in spines
spinetracks = NaN(length(trajectory),1);

% for each spine ROI find overlapping tracks
for k = 1:length(spines)
    xv = spines{k}(:,2); yv = spines{k}(:,1);
    for tj = 1:length(trajectory)
        if trajectory(tj).length >= MTL 
            inout = inpoly2(trajectory(tj).xy*1000, [xv yv]);
            spine = sum(inout)/trajectory(tj).length;   % fraction of track in spine
                       
            % overlapping with spine ROI
            if spine > 0 
                spinetracks(tj) = spine;
            end
        end
    end        
end

% remove spine tracks that were already assigned to dendrite
spinetracks((dendritetracks>0)&(~isnan(spinetracks))) = NaN;

% create structure with only spine associated tracks
spine_trajectory = trajectory(spinetracks>0);

% list of diffusion coefficients in log(um2/s)
Dspine = (log10([spine_trajectory.Dinst]))';

% determine fraction of spine tracks
fraction_spine = length(Dspine)/length(trajectory);

if savefile
    save spine.mat MTL Dspine fraction_spine spinetracks
    save Dspine.txt Dspine -tabs -ascii
end

%% plot results

x = linspace(-4, 1,50);
y1 = hist(Dspine, x); y1=(y1)./sum(y1); 
y2 = hist(Ddendrite, x); y2=(y2)./sum(y2);
% y1 = hist(Dspine, x); 
% y2 = hist(Ddendrite, x); 
figure('Color', 'white');hold on;
plot(x, y1, 'Color',[0.3 0.5 0.75],'LineWidth',2); plot(x, y2, 'Color',[0.9 0.5 0.3], 'LineWidth', 2);
y1=y1.'; y2=y2.';
set(gca,'XLim',[-4 1]);xlabel('log D (um^2/s)'); ylabel('Relative frequency (fractions)');
legend('Spine','Dendrite');
set(gcf,'Position', [100 200 500 400]); 
if savefile
    export_fig 'Hist_spine_dendrite.png'
end

figure('Color', 'white');
c = categorical({'spine','dendrite'});
c = reordercats(c,{'spine' 'dendrite'});
y = [fraction_spine fraction_dendrite];
bar(c,y);ylabel('Fraction of trajectories')
set(gcf,'Position', [650 200 250 400]); 
if savefile
    export_fig 'Fraction_spine_dendrite.png'
end

spine = find(spinetracks > 0);
dendrite = find(dendritetracks > 0);
Dinst = [trajectory.Dinst].';
Dinst = log10(Dinst);

group = zeros(length(spinetracks),1); 
group(spine) = 1; group(dendrite) = 2;
figure('Color', 'white'); boxplot(Dinst, group, 'Labels', {'(not assigned)', 'spine', 'dendrite'})
set(gcf,'Position', [950 200 250 400]); 
ylabel('log D (um^2/s)');
if savefile
    export_fig 'Boxplot_spine_dendrite.png'
end