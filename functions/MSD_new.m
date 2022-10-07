function [M] = MSD_new(xyf,delta) 

% calculate MSD for trajectory
% input:  [x y framenumbers] , dt
% output: [MSD stdev n timelag] for every timelag
%         starting MSD at point 0,0
% compatible with tracks missing frames (e.g. blinking) as it is based
% on frame numbers
%
% modified by Manon Westra, Aug 2021

% Inspired by
% https://nl.mathworks.com/matlabcentral/answers/82268-perform-action-on-small-subset-of-matrices
% https://stackoverflow.com/questions/7489048/calculating-mean-squared-displacement-msd-with-matlab
% &
% "APM_GUI: analyzing particle movement on the cell membrane and determining confinement"
% S.A. Menchon, M.G. Martin and C.G. Dotti

ad = size(xyf);
squaredDisplacement = cell(xyf(ad(1),3)-xyf(1,3)+2,1);

for ai = 1:xyf(ad(1),3)-xyf(1,3)+2
    v(ai) = 0;
end

for an = 1:ad(1)-1
    for ai = 1:ad(1)-an
        dt = xyf(ai+an,3)-xyf(ai,3);
        squaredDisplacement{dt+1} = [squaredDisplacement{dt+1} ;((xyf(ai+an,1)-xyf(ai,1))^2 + (xyf(ai+an,2)-xyf(ai,2))^2)];
        v(dt+1) = v(dt+1)+1;
    end
end

at = 0:delta:(xyf(end,3)-xyf(1,3)+1)*delta;
ann = 1;

meansquaredDisplacement(1,:) = [0 0 0 0]; % MSD curve goes through 0,0

for ai = 1:xyf(ad(1),3)-xyf(1,3)+1
    if v(ai+1)>0
        ann = ann+1;
        meansquaredDisplacement(ann,1) = mean(squaredDisplacement{ai+1});   % mean
        meansquaredDisplacement(ann,2) = std(squaredDisplacement{ai+1});    % stdev
        meansquaredDisplacement(ann,3) = length(squaredDisplacement{ai+1}); % n
        meansquaredDisplacement(ann,4) = at(ai+1);   % time lags
    end
end

M = meansquaredDisplacement;

