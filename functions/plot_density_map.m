function im = plot_density_map(x, y, pxsize, renderpx, show)
%
% from Harold?
% modified by Manon to fix colormap and remove saving

    xx = ceil(x .* (pxsize/renderpx))+1;
    yy = ceil(y .* (pxsize/renderpx))+1;

    im = zeros(max(yy), max(xx));

    for n = 1:length(xx)
        im(yy(n), xx(n)) = im(yy(n), xx(n)) + 1;
    end
   
    
    % map = flipud(rot90(im));
    if show ~ 0;
        density_map = figure;  
        imshow(im,[1 15]); colormap(gca,hot);
    end

