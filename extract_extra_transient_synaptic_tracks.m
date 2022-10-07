%% extract synaptic and extrasynaptic tracks
% based on synaptic mask made with 'make_psd_mask.m'
% extract_confined is needed to define Din and Dout, otherwise comment that
% part out
% INPUT:
%   trajectories.mat
%   PSDmask.mat
% OUTPUT:
%   synaptic.mat
%       diffusion coefficients of the different groups
%       dwell time index
%       fraction synaptic / extrasynaptic
%   updated trajectory struct with synaptic, Dsynaptic, Dextrasynaptic, 
%       inout, entries, exits, Din, Dout, PSD per track
%   histogram of diffusion coefficients for different groups
%
% Manon Westra, 2022


clear;clc;
load trajectories.mat
load PSDmask.mat

extract_extra_transient_synaptic_tracks_param  % load parameters

trajectory_all = trajectory;                %keep original trajectory struct
trajectory_length_all = trajectory_length;  %keep original trajectory_length

% To use only mobile tracks for this analysis
% trajectory = trajectory(mobile);
% trajectory_length = trajectory_length(mobile);

% clear old values
for tj = 1:length(trajectory)
        trajectory(tj).synaptic = [];
        trajectory(tj).Dsynaptic = [];
        trajectory(tj).Dextrasynaptic = [];

        trajectory(tj).inout = [];
        trajectory(tj).entries = [];
        trajectory(tj).exits = [];
        trajectory(tj).Din = [];
        trajectory(tj).Dout = [];
        trajectory(tj).PSD = [];
end
        
% for each PSD find overlapping tracks
for k = 1:length(B)
    xv = B{k}(:,2); yv = B{k}(:,1);
    for tj = 1:length(trajectory)
        if trajectory(tj).length >= MTL 
            x = trajectory(tj).xy(:, 1).*1000; y = trajectory(tj).xy(:, 2).*1000; % convert to nm
            inout = inpolygon (x, y, xv, yv);
            synaptic = sum(inout)/trajectory(tj).length;
                       
            % extrasynaptic but some overlap with PSD
            if synaptic > 0 && synaptic < MFS
                trajectory(tj).synaptic = sum(inout)/trajectory(tj).length;   % fraction of time in synapse
                Dextrasynaptic = trajectory(tj).Dinst;
                trajectory(tj).Dextrasynaptic = Dextrasynaptic;
                trajectory(tj).PSD = k;
                trajectory(tj).inout = inout;
                % determine # entries, # exits
                trans = diff(inout);
                entries = numel(find(trans == 1));
                exits = numel(find(trans == -1));
                trajectory(tj).entries = entries;
                trajectory(tj).exits = exits;
                % determine diffusion in and outside synapse, only possible after running confinement script
                Dr = [trajectory(tj).Dr].';
                Din = mean(Dr(inout == 1));
                Dout = mean(Dr(inout == 0));
                trajectory(tj).Din = Din;
                trajectory(tj).Dout = Dout;
            end
            
            % synaptic
            if synaptic >= MFS
                trajectory(tj).synaptic = sum(inout)/trajectory(tj).length;   % fraction of time in synapse
                Dsynaptic = trajectory(tj).Dinst;
                trajectory(tj).Dsynaptic = Dsynaptic;
                trajectory(tj).PSD = k;
                trajectory(tj).inout = inout;
                                % determine # entries, # exits
                trans = diff(inout);
                entries = numel(find(trans == 1));
                exits = numel(find(trans == -1));
                trajectory(tj).entries = entries;
                trajectory(tj).exits = exits;
                % determine diffusion in and outside synapse, only possible after running confinement script
                Dr = [trajectory(tj).Dr].';
                Din = mean(Dr(inout == 1));
                Dout = mean(Dr(inout == 0));
                trajectory(tj).Din = Din;
                trajectory(tj).Dout = Dout;
            end           
        end
    end        
end

%leftovers with same minimal tracks lengths are extrasynaptic without overlap with PSD:
for tj = 1:length(trajectory)
            if  trajectory(tj).length >= MTL && isempty (trajectory(tj).synaptic)
                Dextrasynaptic = trajectory(tj).Dinst;
                trajectory(tj).synaptic = 0;
                trajectory(tj).Dextrasynaptic = Dextrasynaptic;
                Dr = [trajectory(tj).Dr].';
                trajectory(tj).Dout = mean(Dr);
            end
end

% average fraction of time spent in synapse
fraction_synaptic = [trajectory.synaptic].';
dw = fraction_synaptic > 0;
dwell_time_index = fraction_synaptic(dw);
mean_dti = mean(dwell_time_index);

% list of diffusion coefficients in log(um2/s)
Dsynaptic = (log10([trajectory.Dsynaptic]))';
Dextrasynaptic = (log10([trajectory.Dextrasynaptic]))';
% Dtransient = (log10([trajectory.Dtransient]))';

Din = (log10([trajectory.Din]))';
Dout = (log10([trajectory.Dout]))';

x = linspace(-5, 1,100);
y1 = hist(Dsynaptic, x); y1=(y1)./sum(y1); 
% y2 = hist(Dtransient, x); y2=(y2)./sum(y2);
y3 = hist(Dextrasynaptic, x); y3=(y3)./sum(y3);
y4 = hist(Din, x); y4=(y4)./sum(y4);
y5 = hist(Dout, x); y5=(y5)./sum(y5);
% y1 = hist(Dsynaptic, x); 
% y2 = hist(Dextrasynaptic, x);
% y3 = hist(Dtransient, x); 
figure('Color', 'white');hold on;
plot(x, y1, 'red'); plot(x, y3, 'green');plot(x, y4, 'cyan'); plot(x, y5, 'blue');
y1=y1.';  y3 = y3.';
set(gca,'XLim',[-5 2]);xlabel('log D [um^2/s]');
legend('Synaptic', 'Extrasynaptic', 'Din', 'Dout');
if savefile
    export_fig 'Hist_syn_extra.png'
    save Dsynaptic.txt Dsynaptic -tabs -ascii
    save Dextrasynaptic.txt Dextrasynaptic -tabs -ascii
%     save Dtransient.txt Dtransient -tabs -ascii
end

% Determine fraction synaptic, transient and extrasynaptic
fraction_syn = length(Dsynaptic)/length(trajectory);
fraction_extra = length(Dextrasynaptic)/length(trajectory);
% fraction_trans = length(Dtransient)/length(trajectory);

% fraction_syn_syntrans = length(Dsynaptic)/(length(Dsynaptic)+length(Dtransient));
% fraction_trans_syntrans = length(Dtransient)/(length(Dsynaptic)+length(Dtransient));

% figure('Color', 'white');
% c = categorical({'Synaptic','Extrasynaptic'});
% c = reordercats(c,{'Synaptic' 'Extrasynaptic'});
% y = [fraction_syn fraction_extra];
% bar(c,y);
% if savefile
%     export_fig 'Fraction_syn_extra.png'
% end

% syn = find(fraction_synaptic >= MFS);
% % trans = find((fraction_synaptic > 0) & (fraction_synaptic < MFS));
% extra = find(fraction_synaptic < MFS);
% Dinst = [trajectory.Dinst].';
% Dinst = log10(Dinst);

% group = zeros(length(fraction_synaptic),1); 
% group(syn) = 1; group(extra) = 3;
% figure('Color', 'white'); boxplot(Dinst, group); ylabel('log D [um^2/s]');
% if savefile
%     export_fig 'Boxplot_syn_extra.png'
% end

if savefile
    save('trajectories.mat', 'trajectory', 'trajectory_all', 'trajectory_length_all', '-append')
    save synaptic.mat MTL MFS Dsynaptic Dextrasynaptic fraction_syn fraction_extra ...
        dwell_time_index Din Dout mean_dti
end

%% create structure with only synapse associated tracks
synaptic = [trajectory.synaptic].';
synaptic_tracks = find(synaptic ~= 0);
tracks = trajectory(synaptic_tracks);
