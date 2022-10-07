%% fix when working with older trajectory.mat files that are in nm instead of um

load trajectories.mat trajectory cntr_immobile

%change nm to um
for tj = 1:length(trajectory)
    trajectory(tj).xy = trajectory(tj).xy ./1000;
    trajectory(tj).msd = trajectory(tj).msd ./1000000;
    trajectory(tj).Dinst = trajectory(tj).Dinst ./1000000;    
end

cntr_immobile(:,2:3) = cntr_immobile(:,2:3) ./ 1000;

save('trajectories.mat','trajectory','cntr_immobile', '-append')
