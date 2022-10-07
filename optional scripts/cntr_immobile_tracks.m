%% Retrieve center coordinates of immobile tracks
% for older analyses where this was not yet included in
% convertResults2Trjectory script

load trajectories.mat trajectory immobile

tracks = trajectory(immobile);
cntr_immobile = NaN(length(tracks),3);

for tj = 1:length(tracks)
    xy = tracks(tj).xy;
    mx = mean(xy(:, 1));
    my = mean(xy(:, 2));
    k = immobile(tj);
    cntr_immobile(tj,:) = [k mx my];
end

save('trajectories.mat','cntr_immobile', '-append')