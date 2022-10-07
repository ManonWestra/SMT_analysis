function val = gaussC(x, y, sigma, center)

% from MATLAB Answers, changed line 5 (2*sigma) --> (2*sigma^2)
% https://nl.mathworks.com/matlabcentral/answers/13020-2d-gaussian-function

xc = center(1);
yc = center(2);
exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma^2);
val       = (exp(-exponent));    