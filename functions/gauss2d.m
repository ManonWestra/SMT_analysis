function mat = gauss2d(mat, sigma, center)

% copied from Matlab Answers
% https://nl.mathworks.com/matlabcentral/answers/13020-2d-gaussian-function

gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, sigma, center);