function circle_color(x,y,r,col)
% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% 0.01 is the angle step, bigger values will draw the circle faster but
% you might notice imperfections (not very smooth)
% copied from https://nl.mathworks.com/matlabcentral/answers/3058-plotting-circles

% March 2018 modified by Manon to assign a color to the circle

ang = 0:0.01:2*pi; 
xp = r*cos(ang);
yp = r*sin(ang);
plot(x+xp,y+yp, 'color', col);
end