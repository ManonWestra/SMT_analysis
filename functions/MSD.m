function M = MSD(xy) 

% calculate MSD for trajectory
% input: xy - 

% Inspired by?
% https://nl.mathworks.com/matlabcentral/answers/82268-perform-action-on-small-subset-of-matrices
% https://stackoverflow.com/questions/7489048/calculating-mean-squared-displacement-msd-with-matlab

for dt = 1:length(xy) -1 
    deltaCoords = xy(1+dt:end,1:2) - xy(1:end-dt,1:2);
    squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2
    msd(dt) = mean(squaredDisplacement); %# average
end

M = [0 msd];
