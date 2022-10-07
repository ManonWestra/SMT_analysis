%% Segments PSDs imaged using near-TIRF widefield imaging
% Manon updated Aug 2019 to make the results similar for images from
% different microscopes (different pixel sizes) making use of subsampling

%Update 20200326 Jelmer
    %Easy option to fill in lower and upper limits for widefield display

%Update Dec 2019
    %Saves used Parameters nsub, pkfnd_thresh and thresh_param  
    %Saves PSDpeaks as detected in part one (PSDPeak_xy)
    %Saves all data in one .mat file

clear; close;

make_psd_mask_param  % load parameters

% RESAMPLING PSD IMAGE
% create or load max projection image
if maxroi == 0
    imagelist = dir(filename);
    fileinfo = imfinfo(imagelist.name);
    imsizerow = fileinfo(1).Height; 
    imsizecol = fileinfo(1).Width;
    imsize = [imsizerow imsizecol];
    frames = size(fileinfo,1);
    PSDRaw = MaxProjection(imagelist.name, imsize, frames);
elseif maxroi == 1
    imagelist = dir(filename);
    PSDRaw = imread(imagelist.name);
end

% interpolate
PSDRawsub = imresize(PSDRaw,nsub);

% PSD SEGMENTATION
imsizerow = size(PSDRawsub,1); 
imsizecol = size(PSDRawsub,2);

halfpixel = pixel/2;
pixelsz = pixel/nsub;   % pixel size after resampling

% input bpass function to filter raw image
bpass_lnoise = 2*nsub;
bpass_lobject = round(4*nsub);

% input pkfnd function to detect peaks in PSD image
pkfnd_size = round(250/pixelsz);   % play with this, diameter of blob divided by pixel size to get in pixels
pkfnd_thresh = input('define threshold for peaks: ');

cntrd_region_input = 700/pixelsz;                      % crop region size, length of square in subsampled pixels 
cntrd_region = 2*floor((cntrd_region_input)/2)+1;      % crop region size, should be odd integer for indexing

% Find centers of PSDs 
edge = 1;
PSDFiltered = bpass(PSDRawsub, bpass_lnoise, bpass_lobject);
PSDPeak = pkfnd(PSDFiltered,pkfnd_thresh, pkfnd_size);

[npks,junk]=size(PSDPeak);
    for i=npks:-1:1  % check if too near the border
        if  PSDPeak(i,1)<cntrd_region*edge||PSDPeak(i,1)>imsizecol-cntrd_region*edge || ...
            PSDPeak(i,2)<cntrd_region*edge||PSDPeak(i,2)>imsizerow-cntrd_region*edge 
            PSDPeak(i,:)=[];
        end
    end

% change PSDPeak values into x y coordinates (nm) instead of row column
x = linspace(-halfpixel+0.5*pixelsz,size(PSDRaw,2)*pixel-halfpixel-0.5*pixelsz,size(PSDRaw,2)*nsub);
y = linspace(-halfpixel+0.5*pixelsz,size(PSDRaw,1)*pixel-halfpixel-0.5*pixelsz,size(PSDRaw,1)*nsub);
xpeaks = x(PSDPeak(:,1))';
ypeaks = y(PSDPeak(:,2))';
PSDPeak_xy = [xpeaks ypeaks];

fig = figure('name','TIRF raw image with detected PSDs');
% correct origin image and axes image to nm scale
RI = imref2d(size(PSDRawsub));
RI.XWorldLimits = [(0-halfpixel) (size(PSDRaw,2)*pixel-halfpixel)];
RI.YWorldLimits = [(0-halfpixel) (size(PSDRaw,1)*pixel-halfpixel)];

% plot PSD peak centers + numbers on raw PSD image
imshow(PSDRawsub,RI,[lower_limit_image upper_limit_image]);axis image; hold on;

% or plot on PSDFiltered, where peaks are based on
% imshow(PSDFiltered,RI,[0 200]);axis image; hold on;
 
scatter(PSDPeak_xy(:,1), PSDPeak_xy(:,2),10, [0 0 0], 'filled');
for i = 1:length(PSDPeak)
        text(PSDPeak_xy(i,1)+(5*nsub),PSDPeak_xy(i,2), num2str(i),'Color','r',...
            'FontSize',10);
end

if savefile
    figname = fullfile(cd,[imagelist.name,'_detectedPSDs.fig']);
    saveas(fig, figname);
end   

disp('PSD segmentation done and saved');
%% Individualized PSD filtering

thresh_param = 0.50;   % thresholding parameter, was set as 0.5, to be FWHM-like

% Reconstruct a mask containing the thresheld PSDs
fig = figure('name','PSD mask');
PSDMask = repmat(uint16(0), [imsizerow imsizecol]);
for pkcnt = 1:length(PSDPeak)
    regim=PSDRawsub(PSDPeak(pkcnt,2)-(cntrd_region-1)/2:PSDPeak(pkcnt,2)+(cntrd_region-1)/2,PSDPeak(pkcnt,1)-(cntrd_region-1)/2:PSDPeak(pkcnt,1)+(cntrd_region-1)/2);
%     minint = max([regim(1,1) regim(1,end) regim(end, 1) regim(end,end)]); %corner with max intensity is designated min
    regimlin=double(reshape(regim,1,[]));
    bk=sort(regimlin);
    minint = mean(bk(1:2*cntrd_region));        % mean of (2*cntrd_region) values is designated min 
%     maxint = regim((cntrd_region-1)/2, (cntrd_region-1)/2);
    maxint = mean(bk(end-5:end));               % mean of 6 highest values is designated max
    PSD_thres = minint + (maxint-minint).* thresh_param;       
    PSDMask(PSDPeak(pkcnt,2)-(cntrd_region-1)/2:PSDPeak(pkcnt,2)+(cntrd_region-1)/2,PSDPeak(pkcnt,1)-(cntrd_region-1)/2:PSDPeak(pkcnt,1)+(cntrd_region-1)/2) = ...
        uint16(regim>=PSD_thres);
end

% FURTHER FILTER PSDS BASED ON AREA
% count and label the number of connected components in binary image (i.e.PSDs)
[bw2,n]= bwlabeln(PSDMask);
%figure; imagesc(bw2); axis image;

% determining area
% fill any holes, so that regionprops can be used to estimate
% the area enclosed byhold on each of the boundaries
bw2 = imfill(bw2,'holes');
S = regionprops(bw2, 'Area');

% remove small objects 
P = 0.02*1e6/pixelsz^2; %  WAS 2, BASED ON:  .0188 um^2, minimum PSD area found on EM, ~ 1.875 100nm-pixels
bw2 = ismember(bw2, find([S.Area] > P));

% Use real filtered image as mask, where boundaries are also based on
imshow(bw2,RI,[0 1]);axis image;
if savefile
    figname = fullfile(cd,[imagelist.name,'_PSDMask.fig']);
    saveas(fig, figname);
end

load GreyColormap_binary
colormap(gca,mymap);
if savefile
    figname = fullfile(cd,[imagelist.name,'_PSDMask_grey.fig']);
    saveas(fig, figname);
end

% store centroids and major axis
s  = regionprops(bw2, 'centroid','MajorAxisLength');
%centroids = cat(1, s.Centroid); save 'centroids.mat' centroids;  % WARNING: not correctly adjusted because you need integers for indexing
diameters = cat(1, s.MajorAxisLength); diameters = diameters*pixelsz; 

% find the boundaries
[B,L] = bwboundaries(bw2,'noholes');

% change B to x y coordinates instead of row column
x = linspace(-halfpixel+0.5*pixelsz,size(PSDRaw,2)*pixel-halfpixel-0.5*pixelsz,size(PSDRaw,2)*nsub);
y = linspace(-halfpixel+0.5*pixelsz,size(PSDRaw,1)*pixel-halfpixel-0.5*pixelsz,size(PSDRaw,1)*nsub);

% show boundaries of detected PSDs
fig = figure('name','PSD mask: final stage of filtering');
imshow(PSDRawsub,RI,[lower_limit_image upper_limit_image]);axis image; hold on;

% plot the boundaries
for k = 1:length(B)
  boundary = B{k};
  xboundary = x(boundary(:,2))';
  yboundary = y(boundary(:,1))';
  boundary = [yboundary xboundary];
  B{k} = boundary;
  plot(boundary(:,2), boundary(:,1), 'LineWidth', 1, 'Color', 'red');
  text(boundary(1,2)+(5*nsub),boundary(1,1), num2str(k),'Color','b',...
       'FontSize',10);
end
axis image

if savefile
    figname = fullfile(cd,[imagelist.name,'_PSDboundaries.fig']);
    saveas(fig, figname);
end

% Calculate area of detected PSDs and save stats
stats = regionprops(L,'Area','MajorAxisLength','MinorAxisLength','Eccentricity','Perimeter');
areaPSDsub = [stats.Area]';
areaPSD = areaPSDsub*pixelsz^2/1e6; 
areaMean = mean(areaPSD);
save PSDmask.mat B stats nsub pixel pkfnd_thresh thresh_param PSDPeak_xy areaPSD areaMean diameters

disp('Individualized PSD filtering done and saved');