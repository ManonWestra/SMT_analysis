function [data] = get_roi_coords(locs)
% get coordinates in ROI together with their frame number
%
% from Harold?
% modified by Manon to save also frame number with xy coordinates

figure('Color', 'white');
scatter(locs(:,1), locs(:,2),1,'red','filled'); 
axis equal;
h = imrect;
p = getPosition(h);

a = find(locs(:,1) > p(1) & locs(:,1) < p(1) + p(3));
b = find(locs(a,2) > p(2) & locs(a,2) < p(2) + p(4));
x = locs(a(b),1);
y = locs(a(b),2);
frame = locs(a(b),3);
data = [x y frame];
mean(x)
mean(y)

